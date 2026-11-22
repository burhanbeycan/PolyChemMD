"""Exercise the local application through Chromium, using a fresh job store.

Install Studio plus playwright, then run 'python -m playwright install chromium'.
An existing Chromium binary can be selected with --chromium. Use --with-md
after installing LAMMPS, MDAnalysis and the reference IFF case.
"""

import argparse
import json
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import time
from urllib.request import urlopen


def read_json(url):
    with urlopen(url, timeout=3) as response:
        return json.load(response)


def check(page, url, directory, with_md):
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(url)
    page.get_by_role("heading", name="Sistem oluştur").wait_for()

    def submit(button, expected="succeeded"):
        def matches(response):
            return (
                response.url == url + "/api/jobs" and response.request.method == "POST"
            )

        with page.expect_response(matches) as response:
            button.click()
        job_id = response.value.json()["id"]
        deadline = time.monotonic() + 120
        while time.monotonic() < deadline:
            job = read_json(f"{url}/api/jobs/{job_id}")
            if job["status"] in {"succeeded", "failed"}:
                break
            page.wait_for_timeout(250)
        if job["status"] != expected:
            raise AssertionError(f"{job['kind']} {job['status']}: {job['error']}")
        page.locator("#refresh-jobs").click()
        status = "Tamamlandı" if expected == "succeeded" else "Başarısız"
        page.locator("#job-detail").get_by_text(status, exact=True).wait_for()
        return job

    job = submit(page.locator("button[value=packing]"))
    with page.expect_download() as download:
        page.locator('a[download="packed.xyz"]').click()
    xyz = directory / "coordinates.xyz"
    download.value.save_as(xyz)
    with xyz.open() as stream:
        assert int(stream.readline()) == job["result"]["atom_count"]
    page.get_by_text("3B yapıyı aç", exact=True).click()
    page.wait_for_function(
        "document.querySelector('#view-source').textContent.includes('Paketleme')"
    )
    canvas = page.locator("#scene")
    initial = canvas.evaluate("c => c.toDataURL()")
    canvas.press("ArrowRight")
    page.wait_for_function(
        "initial => document.querySelector('#scene').toDataURL() !== initial",
        arg=initial,
    )

    page.get_by_role("button", name="Deneyle karşılaştır").click()
    xrd = submit(page.locator("#reference-run"))
    assert abs(xrd["result"]["summary"]["mean"] + 0.172) < 1e-9
    page.get_by_role("button", name="Atlas kayıtları").click()
    page.locator("#record-search").fill("Xing")
    page.wait_for_function(
        "document.querySelector('#record-count').textContent.startsWith('2 /')"
    )
    page.get_by_role("button", name="Hesaplama yöntemleri").click()
    page.locator("#plan-form button").click()
    page.locator("#plan-result .plan-job").nth(1).wait_for()

    if with_md:
        page.locator("#md-form [name=md_steps]").fill("4")
        md = submit(page.locator('#md-form button[type="submit"]'))
        assert md["result"]["last_sampled_step"] == 4
        trajectory = submit(page.get_by_role("button", name="Yörüngeyi analiz et"))
        assert trajectory["result"]["n_frames"] == md["result"]["sampled_frames"]
        assert 0 <= trajectory["result"]["mean_water_O_contact_fraction"] <= 1

    page.get_by_role("button", name="Deneyle karşılaştır").click()
    page.locator("#experiment-form [name=assay]").select_option("optical")
    for key, value in {
        "batch": "QA1",
        "medium": "water",
        "polymer": "PCE",
        "A0": "0.1",
        "A24": "0.05",
        "B0": "0.15",
        "B24": "0.01",
    }.items():
        page.locator(f"#assay-grid [data-key={key}]").fill(value)
    rejected = submit(page.locator('#experiment-form button[type="submit"]'), "failed")
    assert "baseline" in rejected["error"]
    page.set_viewport_size({"width": 390, "height": 844})
    page.get_by_role("button", name="Sistem oluştur").click()
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    assert not errors, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromium", help="Path to an installed Chromium executable")
    parser.add_argument("--with-md", action="store_true")
    args = parser.parse_args()
    from playwright.sync_api import sync_playwright

    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        url = f"http://127.0.0.1:{port}"
        cmd = [
            sys.executable,
            "-m",
            "polychemmd.studio",
            "--data-dir",
            tmp,
            "serve",
            "--port",
            str(port),
        ]
        with subprocess.Popen(cmd) as server:
            try:
                for _ in range(100):
                    if server.poll() is not None:
                        raise RuntimeError("Studio server exited during startup")
                    try:
                        read_json(url + "/api/status")
                        break
                    except OSError:
                        time.sleep(0.15)
                else:
                    raise RuntimeError("Studio server did not become ready")
                with sync_playwright() as playwright:
                    browser = playwright.chromium.launch(executable_path=args.chromium)
                    page = browser.new_page(
                        viewport={"width": 1480, "height": 950}, accept_downloads=True
                    )
                    check(page, url, directory, args.with_md)
                    browser.close()
            finally:
                server.terminate()
    print("Chromium Studio workflow passed.")


if __name__ == "__main__":
    main()
