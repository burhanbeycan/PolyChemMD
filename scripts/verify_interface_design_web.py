"""Exercise controlled microstructures, packing, topology generation and packing comparisons.

Run with Python from the Studio environment and Playwright, using --chromium
to select a local browser. --with-md adds generated topology and dynamics to
the checks and requires pinned IFF source data, Moltemplate and LAMMPS.
"""

import time

from playwright.sync_api import expect

from verify_studio_web import main, read_json


def verify(page, url, directory, with_md, *, screenshot=None):
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(url)
    page.get_by_role("heading", name="Sistem oluştur").wait_for()
    detail = page.locator("#job-detail")

    def run(button):
        with page.expect_response(
            lambda response: (
                response.url == url + "/api/jobs" and response.request.method == "POST"
            )
        ) as submitted:
            button.click()
        job_id = submitted.value.json()["id"]
        deadline = time.monotonic() + 240
        while time.monotonic() < deadline:
            job = read_json(f"{url}/api/jobs/{job_id}")
            if job["status"] not in ("running", "queued"):
                break
            page.wait_for_timeout(250)
        if job["status"] != "succeeded":
            raise RuntimeError(f"{job['kind']} {job['status']}: {job['error']}")
        page.locator("#refresh-jobs").click()
        detail.get_by_text("Tamamlandı", exact=True).wait_for()
        return job

    page.locator("#designed-pva").click()
    builder = page.locator("#system-form")
    assert builder.locator('[name="salt_M"]').input_value() == "0"
    builder.locator('[name="acetate_pattern"]').select_option("explicit")
    positions = builder.locator('[name="acetate_positions"]')
    positions.fill("2, 2, 9")
    assert not positions.evaluate("e => e.validity.valid")
    builder.locator("[name=surface_correlation_top]").fill("0.4")
    builder.locator("[name=surface_correlation_bottom]").fill("-0.2")
    positions.fill("3, 4, 9")
    assert positions.evaluate("e => e.validity.valid")
    packed = run(builder.locator('button[value="packing"]'))
    detail.get_by_role("heading", name="Gerçekleşen asetat dizisi").wait_for()
    assert detail.locator(".sequence-result .acetate").all_text_contents() == [
        "3",
        "4",
        "9",
    ]
    assert (
        detail.get_by_role("table", name="Yüzey komşuluğu").locator("tbody tr").count()
        == 2
    )

    detail.get_by_role("button", name="3B yapıyı aç").click()
    page.locator(".face-map").first.wait_for()
    assert page.locator(".face-map").count() == 2
    assert page.locator(".face-map circle").count() == 96
    page.locator("#toggle-water").click()
    expect(page.locator("#toggle-water")).to_have_attribute("aria-pressed", "false")
    if screenshot:
        page.locator("#system-form details").first.evaluate("node => node.open = true")
        page.locator("#notice").wait_for(state="hidden", timeout=10000)
        page.evaluate("window.scrollTo(0,0)")
        page.screenshot(path=screenshot, animations="disabled")
    page.get_by_role("button", name="Hesaplama kuyruğu").click()
    if with_md:
        detail.get_by_role("button", name="Bu yapıdan topoloji üret").click()
        generator = page.locator("#topology-form")
        expect(generator.locator('[name="source_job"]')).to_have_value(packed["id"])
        generated = run(generator.get_by_role("button", name="Topolojiyi üret"))
        assert generated["result"]["atoms"] == packed["result"]["atom_count"]
        assert (
            generated["result"]["charge_e"] == 0
            or abs(generated["result"]["charge_e"]) < 1e-5
        )
        with page.expect_download() as download:
            detail.locator('a[download="typed.data"]').click()
        typed = directory / "built-topology.data"
        download.value.save_as(typed)
        assert typed.stat().st_size > 1000
        detail.get_by_role("button", name="Bu topolojiyle dinamik çalıştır").click()
        registered = page.locator("#registered-form")
        page.wait_for_function(
            "document.querySelector('#registered-form [name=topology_id]').value !== ''"
        )
        assert (
            registered.locator("[name=topology_id]").input_value()
            == generated["result"]["topology_id"]
        )
        registered.locator("[name=steps]").fill("4")
        pilot = run(
            registered.get_by_role("button", name="Parametreli MD işini sırala")
        )
        assert pilot["result"]["steps"] == 4
        assert pilot["result"]["post_minimize_mobile_fmax_kcal_mol_A"] >= 0

    page.get_by_role("button", name="Hesaplama yöntemleri").click()
    bench = page.locator("#benchmark-form")
    for value in ("molpack", "packmol", "confined"):
        bench.locator(f'input[value="{value}"]').uncheck()
    bench.locator("input[name=seeds]").fill("17")
    comparison = run(bench.get_by_role("button"))
    assert len(comparison["result"]["summary"]) == 2
    assert all(row["success_count"] == 1 for row in comparison["result"]["summary"])
    assert (
        detail.get_by_role("table", name="Paketleme sonuçları")
        .locator("tbody tr")
        .count()
        == 2
    )
    page.get_by_role("button", name="Sistem oluştur").click()
    page.get_by_text(
        "Ölçülen Mₙ ve dispersiteden PVA zincirleri oluştur", exact=True
    ).click()
    population = run(page.locator("#ensemble-form button"))
    choice = population["result"]["members"][0]
    detail.get_by_role("button", name="Seçili zinciri kullan").click()
    assert positions.input_value() == ",".join(map(str, choice["acetate_positions"]))
    assert int(builder.locator("[name=dp]").input_value()) == choice["dp"]
    assert positions.evaluate("e => e.validity.valid")
    page.set_viewport_size({"width": 390, "height": 844})
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    assert not errors, errors


if __name__ == "__main__":
    import argparse
    from functools import partial
    from pathlib import Path
    import sys

    captures = argparse.ArgumentParser(add_help=False)
    captures.add_argument("--screenshot", type=Path)
    options, rest = captures.parse_known_args()
    sys.argv[1:] = rest
    main(
        partial(verify, screenshot=options.screenshot),
        __doc__ + "\nUse --screenshot FILE to save the completed interface view.",
    )
