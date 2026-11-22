// Keep recipe data numeric at the API boundary. Blank seed overrides stay absent.
export function integerList(text) {
  if (!text.trim()) return [];
  const tokens = text.split(',').map(value => value.trim());
  return tokens.every(token => /^\d+$/.test(token)) ? tokens.map(Number) : null;
}

function group(form, id, active) {
  const section = form.querySelector('#' + id);
  section.hidden = !active;
  for (const control of section.querySelectorAll('input,select')) control.disabled = !active;
}

export function microstructureFields(form, config) {
  const fields = form.elements;
  for (const key of ['surface_seed', 'sequence_seed', 'tacticity_seed', 'conformer_seed']) {
    if (fields[key].value !== '' && (config.polymer !== 'bare' || key === 'surface_seed')) config[key] = Number(fields[key].value);
  }
  if (fields.packing_seed.value !== '') config.packing_seed = Number(fields.packing_seed.value);
  if (config.pattern === 'correlated') {
    config.surface_correlation_top = Number(fields.surface_correlation_top.value);
    config.surface_correlation_bottom = Number(fields.surface_correlation_bottom.value);
  }
  if (config.polymer !== 'PVA') return;
  if (config.acetate_pattern === 'dyads') {
    const key = fields.dyad_measure.value === 'fraction' ? 'acetate_dyad_fraction' : 'acetate_dyads';
    config[key] = Number(fields[key].value);
  }
  if (config.acetate_pattern === 'explicit') config.acetate_positions = integerList(fields.acetate_positions.value);
}

export function syncMicrostructure(form, config) {
  const fields = form.elements, pva = config.polymer === 'PVA';
  const dyads = pva && fields.acetate_pattern.value === 'dyads';
  const explicit = pva && fields.acetate_pattern.value === 'explicit';
  group(form, 'correlation-config', config.pattern === 'correlated');
  group(form, 'pva-config', pva);
  group(form, 'pce-config', config.polymer === 'PCE');
  group(form, 'dyad-config', dyads);
  group(form, 'dyad-count-label', dyads && fields.dyad_measure.value === 'count');
  group(form, 'dyad-fraction-label', dyads && fields.dyad_measure.value === 'fraction');
  group(form, 'acetate-positions-label', explicit);
  fields.acetate_dyads.required = dyads;
  fields.acetate_dyad_fraction.required = dyads;
  fields.surface_correlation_top.required = fields.surface_correlation_bottom.required = config.pattern === 'correlated';
  const n = config.dp;
  const k = Math.floor((1 - Number(fields.hydrolysis.value)) * n + 0.5 + 1e-12);
  const minimum = Math.max(0, 2*k - n - 1), maximum = Math.max(0, k - 1);
  fields.acetate_dyads.min = String(minimum);
  fields.acetate_dyads.max = String(maximum);
  const positions = integerList(fields.acetate_positions.value);
  const valid = positions !== null && positions.length === k && positions.every(v => v >= 1 && v <= n) && new Set(positions).size === k;
  fields.acetate_positions.setCustomValidity(explicit && !valid ? 'Hidroliz kesrine uygun ' + k + ' farklı tekrar konumunu (1–' + n + ') girin.' : '');
  form.querySelector('#sequence-budget').textContent = Number.isInteger(n) && n > 0
    ? n + ' tekrarda ' + k + ' asetat. Mümkün komşu OAc–OAc sayısı ' + minimum + '–' + maximum + '.' : '';
}

function cell(tag, value) {
  const node = document.createElement(tag);
  node.textContent = value;
  return node;
}

function appendTable(parent, headers, rows, label) {
  if (!rows.length) return;
  const frame = document.createElement('div');
  frame.className = 'table-scroll';
  const table = document.createElement('table');
  table.setAttribute('aria-label', label);
  const head = table.createTHead().insertRow(), body = table.createTBody();
  for (const text of headers) head.append(cell('th', text));
  for (const values of rows) {
    const row = body.insertRow();
    for (const value of values) row.append(cell('td', value));
  }
  frame.append(table); parent.append(frame);
}

export function drawStructureStats(parent, shape, fmt) {
  const polymer = shape.polymer, stats = polymer?.sequence_statistics;
  if (stats && Array.isArray(polymer.acetate_residues_1based)) {
    parent.append(cell('h3', 'Gerçekleşen asetat dizisi'));
    const list = document.createElement('ol'), positions = new Set(polymer.acetate_residues_1based);
    list.className = 'sequence-result';
    for (let i = 1; i <= polymer.degree_of_polymerization; i++) {
      const acetate = positions.has(i);
      const item = cell('li', String(i));
      item.title = item.ariaLabel = i + '. tekrar · ' + (acetate ? 'asetat' : 'alkol');
      if (acetate) item.className = 'acetate';
      list.append(item);
    }
    parent.append(list, cell('p', 'Asetat konumları: ' + (Array.from(positions).join(', ') || 'yok') +
      ' · Komşu OAc–OAc: ' + fmt(stats.acetate_dyads) + ' · Çift kesri: ' + fmt(stats.acetate_dyad_fraction)));
  }
  const rows = [];
  for (const [sheet, faces] of (shape.face_occupancy || []).entries()) {
    for (const [face, label] of [['top', 'üst'], ['bottom', 'alt']]) {
      const r = faces[face], counts = r?.neighbor_statistics;
      if (counts) rows.push([String(sheet + 1) + ' / ' + label, fmt(r.n_OH), fmt(r.n_F),
        fmt(counts.requested_correlation), fmt(counts.oh_correlation), fmt(counts.absolute_correlation_error)]);
    }
  }
  if (rows.length) {
    parent.append(cell('h3', 'Yüzey komşuluğu'));
    appendTable(parent, ['Tabaka / yüz', 'OH', 'F', 'Hedef korelasyon', 'Gerçekleşen', '|Fark|'], rows, 'Yüzey komşuluğu');
  }
}

const engineNames = {builtin:'Yerleşik', slit:'Dar galeri', molpack:'MolPack', packmol:'Packmol'};
const caseNames = {aqueous:'Sulu arayüz', confined:'Dar galeri'};
export function drawBenchmark(parent, result, fmt) {
  parent.append(cell('h3', 'Motor karşılaştırması'));
  const rows = (result.summary || []).map(r => [caseNames[r.case] || r.case, engineNames[r.engine] || r.engine,
    fmt(r.success_count), fmt(r.failure_count), fmt(r.unavailable_count), fmt(r.timeout_count), fmt(r.median_success_s), fmt(r.minimum_distance_A)]);
  appendTable(parent, ['Koşul', 'Motor', 'Başarı', 'Hata', 'Kurulu değil', 'Zaman aşımı', 'Medyan süre (s)', 'En kısa temas (Å)'], rows, 'Paketleme sonuçları');
  const checks = (result.rows || []).filter(r => r.reason || r.force_status && r.force_status !== 'not_requested');
  if (checks.length) {
    const d = document.createElement('details'); d.append(cell('summary', 'Koşu bazında kuvvetler ve tanılar'));
    appendTable(d, ['Koşul / motor', 'Tohum', 'Kuvvet durumu', 'İlk Fmax', 'Son Fmax', 'Açıklama'],
      checks.map(r => [caseNames[r.case] + ' / ' + engineNames[r.engine], fmt(r.seed), r.force_status || '—',
        fmt(r.initial_mobile_fmax_kcal_mol_A), fmt(r.final_mobile_fmax_kcal_mol_A), r.reason || r.force_reason || '—']),
      'Hareketli atom kuvvetleri (kcal mol⁻¹ Å⁻¹) ve tanılar');
    d.append(cell('p','Fmax birimi kcal mol⁻¹ Å⁻¹.')); parent.append(d);
  }
}

export function drawFaceMaps(host, atoms, shape, fmt) {
  host.replaceChildren();
  if (!Array.isArray(shape?.face_occupancy)) return false;
  const [lx, ly] = shape.box_A, ns = 'http://www.w3.org/2000/svg';
  for (const [sheet, faces] of shape.face_occupancy.entries()) {
    for (const [face, label] of [['top', 'üst'], ['bottom', 'alt']]) {
      const f = faces[face], stats = f.neighbor_statistics;
      if (!stats) continue;
      const map = document.createElementNS(ns, 'svg');
      map.setAttribute('viewBox', '-1 -1 ' + (lx+2) + ' ' + (ly+2));
      map.setAttribute('aria-label', (sheet+1) + '. tabaka ' + label + ' yüzeyinin OH ve F sonlanma konumları');
      map.setAttribute('role', 'img');
      for (const a of atoms) {
        if (a.sheet_index !== sheet || !['O','F'].includes(a.element) ||
            ![face+'_OH',face+'_F'].includes(a.layer)) continue;
        const dot = document.createElementNS(ns, 'circle');
        dot.setAttribute('cx', String(a.x));
        dot.setAttribute('cy', String(ly-a.y));
        dot.setAttribute('r','0.91');
        dot.setAttribute('fill', a.element === 'F' ? '#a96334' : '#197e75');
        map.append(dot);
      }
      const card = document.createElement('div'), desc = document.createElement('div');
      card.className = 'face-map'; desc.append(cell('strong', 'Tabaka '+(sheet+1)+' · '+label+' yüz'),
        cell('small', fmt(f.n_OH)+' OH / '+fmt(f.n_F)+' F'),
        cell('small','Komşuluk ρ = '+fmt(stats.oh_correlation)));
      card.append(map,desc); host.append(card);
    }
  }
  return host.childElementCount > 0;
}
