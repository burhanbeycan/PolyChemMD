import {StructureView} from './viewer.js';
import {drawBenchmark, drawFaceMaps, drawStructureStats, integerList, microstructureFields, syncMicrostructure} from './design.js';
const $=s=>document.querySelector(s);
const $$=s=>Array.from(document.querySelectorAll(s));
const svgNS='http://www.w3.org/2000/svg';
const icons={layers:'M2 7l10-5 10 5-10 5-10-5zm0 5 10 5 10-5M2 17l10 5 10-5',
 folder:'M3 5h6l2 3h10v12H3z',sliders:'M3 6h4m4 0h10M3 12h10m4 0h4M3 18h6m4 0h8M7 3v6m6 0v6m-4 0v6',
 queue:'M3 5h18M3 10h14M3 15h9m4-2 5 4-5 4z',flask:'M9 3h6m-5 0v6l-6 10q-1 2 2 2h12q3 0 2-2l-6-10V3M7 15h10',
 book:'M12 5c-3-2-6-2-10-1v15c4-1 7-1 10 1 3-2 6-2 10-1V4c-4-1-7-1-10 1zm0 0v15',
 terminal:'M3 5h18v14H3z M6 9l3 3-3 3m5 0h5',reset:'M3 11a9 9 0 1 1 2 7M3 5v6h6',
 expand:'M3 9V3h6m6 0h6v6m0 6v6h-6m-6 0H3v-6',arrow:'M3 12h17m-6-6 6 6-6 6',play:'M6 3l14 9-14 9z'};
for(const icon of $$('[data-icon]')){
  const svg=document.createElementNS(svgNS,'svg'),path=document.createElementNS(svgNS,'path');
  svg.setAttribute('viewBox','0 0 24 24');svg.setAttribute('class','icon');svg.setAttribute('aria-hidden','true');
  path.setAttribute('d',icons[icon.dataset.icon]);svg.append(path);icon.append(svg);
}
const element=(tag,text='',cls='')=>{const el=document.createElement(tag);if(text!==undefined)el.textContent=text;if(cls)el.className=cls;return el};
const state={page:'builder',selected:null,jobs:[],jobSignature:'',records:[],atlasLoaded:false,service:false,viewJob:null,
             rules:null,polling:false,jobRequest:0};
const viewer=new StructureView($('#scene'));
const kindLabel={geometry:'Geometri',packing:'Paketleme',md:'Moleküler dinamik',radonpy:'RadonPy',experiment:'Deney analizi',trajectory:'Yörünge analizi',topology:'Topoloji üretimi',packing_benchmark:'Paketleme karşılaştırması',pva_ensemble:'PVA zincir dağılımı'};
const statusLabel={queued:'Sırada',running:'Çalışıyor',succeeded:'Tamamlandı',failed:'Başarısız',partial_failed:'Kısmen başarısız',cancelled:'İptal edildi'};
const f=$('#system-form'), field=name=>f.elements.namedItem(name);
let noticeTimeout;
function notice(message,error=false){const el=$('#notice');el.textContent=message;el.classList.toggle('error',error);el.hidden=false;clearTimeout(noticeTimeout);if(!error)noticeTimeout=setTimeout(()=>{el.hidden=true},5800)}
async function api(path,options){
  const res=await fetch(path,options);
  const json=await res.json().catch(()=>null);
  if(!res.ok)throw Error(typeof json?.detail==='string'?json.detail:('HTTP '+res.status));
  if(json===null)throw Error('Geçerli bir JSON yanıtı alınamadı.');return json;
}
function go(name,push=true){
  if(!['builder','methods','jobs','experiment','atlas'].includes(name))name='builder';
  for(const page of $$('.page'))page.hidden=page.id!=='page-'+name;
  for(const btn of $$('nav button')){if(btn.dataset.page===name)btn.setAttribute('aria-current','page');else btn.removeAttribute('aria-current')}
  state.page=name;if(push&&location.hash!=='#'+name) history.replaceState(null,'','#'+name);
  if(name==='builder')viewer.schedule();
  if(name==='atlas')loadAtlas();
  if(name==='methods'){loadRegistered();loadCampaign()}
  if(name==='jobs')loadJobs();
}
for(const btn of $$('[data-page],[data-go]'))btn.addEventListener('click',()=>go(btn.dataset.page||btn.dataset.go));
addEventListener('hashchange',()=>go(location.hash.slice(1),false));
function systemConfig(){
  const n=k=>Number(field(k).value),v=k=>field(k).value;
  const c={polymer:v('polymer'),geometry:v('geometry'),vacuum_padding_A:n('vacuum_padding_A'),dp:n('dp'),chains:v('polymer')==='bare'?0:n('chains'),
    seed:n('seed'),nx:n('nx'),ny:n('ny'),gap_A:n('width_A'),water_height_A:n('width_A'),oh_top:n('oh_top'),oh_bottom:n('oh_bottom'),
    pattern:v('pattern'),salt_M:n('salt_M'),water_density:n('water_density'),conformer_ff:v('conformer_ff'),packing_engine:v('packing_engine')};
  if(c.polymer==='PVA'){c.hydrolysis=n('hydrolysis');c.acetate_pattern=v('acetate_pattern')}
  if(c.polymer==='PCE'){c.side_dp=n('side_dp');c.graft_every=n('graft_every')}
  microstructureFields(f,c);
  return c;
}
function syncBuilder(){
  const c=systemConfig(),bare=c.polymer==='bare';$('#project-name').textContent='Ti₃C₂Tₓ / '+(bare?'kontrol':c.polymer);
  $('#pva-config').hidden=c.polymer!=='PVA';$('#pce-config').hidden=c.polymer!=='PCE';$('#chain-fields').hidden=bare;
  syncMicrostructure(f,c);
  field('dp').disabled=bare;field('chains').disabled=bare;
  $('#width-label').textContent=c.geometry==='slit'?'Galeri genişliği (Å)':'Su bölgesi (Å)';
  $('#match-name').textContent=c.geometry==='slit'?'Kuru / ıslak bazal tekrar mesafesi':'Yüzeye dik hidrasyon profili';
  $('#match-control').textContent=c.geometry==='slit'?'Aynı durumda polimersiz MXene':'Aynı sonlanmalara sahip çıplak yüzey';
  if(state.viewJob){state.viewJob=null;$('#view-source').textContent='Şematik yerleşim · Ölçeksiz';$('#view-detail').textContent='Şema'}
  $('#site-maps').hidden=true;$('#face-maps').replaceChildren();$('#toggle-water').hidden=true;$('#water-legend').classList.remove('omitted');
  viewer.schematic(c);
}
f.addEventListener('input',syncBuilder);
$('#toggle-water').onclick=e=>{
  viewer.setWaterVisible(!viewer.waterVisible);
  e.currentTarget.setAttribute('aria-pressed',String(viewer.waterVisible));
  $('#water-legend').classList.toggle('omitted',!viewer.waterVisible);
  e.currentTarget.title=viewer.waterVisible?'Suyu gizle':'Suyu göster';viewer.schedule();
};
$('#reset-view').onclick=()=>viewer.reset();
$('#expand-view').onclick=async()=>{
  const card=$('#viewer-card');try{if(document.fullscreenElement)await document.exitFullscreen();else if(card.requestFullscreen)await card.requestFullscreen();else card.classList.toggle('expanded')}
  catch{card.classList.toggle('expanded')} viewer.schedule();
};
addEventListener('keydown',e=>{if(e.key==='Escape')$('#viewer-card').classList.remove('expanded')});
async function safe(fn,trigger){
  if(trigger?.disabled)return;
  if(trigger)trigger.disabled=true;
  try{return await fn()}catch(e){notice(e.message||String(e),true)}finally{if(trigger)trigger.disabled=false}
}
async function submit(kind,config){
  const val=await api('/api/jobs',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({kind,config})});
  const job=val.job||val;state.selected=job.id;state.jobSignature='';go('jobs');await loadJobs(true);
  notice(kindLabel[kind]+' işi kuyruğa alındı.');return job;
}
f.onsubmit=e=>{
  e.preventDefault();
  if(!f.reportValidity())return;
  const kind=e.submitter?.value||'geometry', config=systemConfig();
  if(kind==='geometry')delete config.packing_seed;
  safe(()=>submit(kind,config),e.submitter);
};
function applyDesign(values) {
  for(const [key,value] of Object.entries(values))field(key).value=String(value);
  if('acetate_dyad_fraction' in values)field('dyad_measure').value='fraction';
  syncBuilder();
}
$('#designed-pva').onclick=()=>{
  f.reset();
  applyDesign({polymer:'PVA',geometry:'single',dp:12,hydrolysis:0.75,acetate_pattern:'dyads',acetate_dyads:1,
    pattern:'correlated',surface_correlation_top:0,surface_correlation_bottom:0,salt_M:0,packing_engine:'slit',vacuum_padding_A:12});
  field('acetate_pattern').closest('details').open=true;
  notice('Tuzsuz DP12 / %75 hidroliz örneği hazır. Dizilim ve yüzey hedeflerini değiştirebilirsiniz.');
};
const ensemble=$('#ensemble-form');
ensemble.onsubmit=e=>{
  e.preventDefault();
  const config=Object.fromEntries(Array.from(new FormData(ensemble),([key,value])=>[key,Number(value)]));
  safe(()=>submit('pva_ensemble',config),e.submitter);
};
const benchmark=$('#benchmark-form');
function benchmarkConfig(){
  const pick=name=>Array.from(benchmark.querySelectorAll('[name='+name+']:checked'),node=>node.value);
  return {cases:pick('case'),engines:pick('engine'),seeds:integerList(benchmark.elements.seeds.value),
    evaluate_forces:benchmark.elements.evaluate_forces.checked};
}
function syncBenchmark(){
  const c=benchmarkConfig();
  const valid=c.seeds?.length && c.seeds.length<=8 && c.seeds.every(v=>Number.isSafeInteger(v) && v>=1 && v<=2147483646) &&
    new Set(c.seeds).size===c.seeds.length;
  const msg=!valid?'1–8 farklı pozitif tam sayı tohumunu virgülle ayırın.':
    !c.cases.length || !c.engines.length?'En az bir koşul ve bir motor seçin.':'';
  benchmark.elements.seeds.setCustomValidity(msg);
  $('#benchmark-budget').textContent=msg || (c.cases.length*c.engines.length*c.seeds.length+' paketleme koşusu');
}
benchmark.oninput=syncBenchmark;
benchmark.onsubmit=e=>{
  e.preventDefault();syncBenchmark();
  if(benchmark.reportValidity())safe(()=>submit('packing_benchmark',benchmarkConfig()),e.submitter);
};
syncBenchmark();
const md=$('#md-form');
function mdAvailability(){
  const profile=md.elements.md_profile.value, key=profile==='iff-water'?profile:profile+'-'+md.elements.md_sequence.value;
  $('#md-seq-label').hidden=profile==='iff-water';
  if(!state.iffProfiles || !(key in state.iffProfiles))return;
  $('#md-availability').textContent=state.iffProfiles[key]?'Kaynak dosyaları yerel girdi deposunda bulundu.':'Bu reçete için IFF girdi dosyalarını önce yerel depoya aktarın (kurulum kılavuzu).';
}
md.elements.md_profile.onchange=mdAvailability;md.elements.md_sequence.onchange=mdAvailability;
md.onsubmit=e=>{e.preventDefault();const a=md.elements;safe(()=>submit('md',{md_profile:a.md_profile.value,md_sequence:a.md_sequence.value,md_steps:Number(a.md_steps.value),seed:Number(a.seed.value)}),e.submitter)};
const radon=$('#radon-form');
radon.onsubmit=e=>{e.preventDefault();const a=radon.elements;safe(()=>submit('radonpy',{polymer:a.polymer.value,dp:Number(a.dp.value),seed:Number(a.seed.value)}),e.submitter)};
function toolBadge(ok,version){const txt=ok?(version||'Kurulu'):'Bulunamadı';return element('span',txt,'status-tag'+(ok?'':' warn'))}
async function loadStatus(){
  const shown={rdkit:'RDKit',ase:'ASE',slit:'PolyChemMD · dar galeri',moltemplate:'Moltemplate',molpack:'MolCrafts MolPack',packmol:'Packmol',lammps:'LAMMPS',mdanalysis:'MDAnalysis',radonpy:'RadonPy (bağımsız ortam)',builtin:'PolyChemMD'};
  const statusEl=$('#service');
  try{
    const data=await api('/api/status');state.service=true;state.iffProfiles=data.profiles;mdAvailability();statusEl.lastChild.textContent=' Yerel çalışma alanı';statusEl.classList.remove('offline');
    const el=$('#tools');el.replaceChildren();
    for(const [name,m] of Object.entries(data.tools||{})){
      const row=element('div','','tool-row');row.append(element('span',shown[name]||name),toolBadge(Boolean(m.available),m.version));el.append(row);
    }
  }catch{state.service=false;statusEl.lastChild.textContent=' Servise erişilemiyor';statusEl.classList.add('offline')}
}
function jobStatus(s){return element('span',statusLabel[s]||s,'status-tag'+(s==='failed'?' fail':(['running','queued'].includes(s)?' warn':'')))}
const formatNum=v=>typeof v==='number'?v.toLocaleString('tr-TR',{maximumFractionDigits:6}):String(v??'—');
function showJobs(){
  const el=$('#job-list');el.replaceChildren();
  if(!state.jobs.length){el.append(element('p','Henüz iş yok. Sistem oluştur ekranından ilk geometriyi hazırlayabilirsiniz.','empty'));return}
  for(const j of state.jobs){
    const b=element('button','','job-card'+(j.id===state.selected?' selected':''));b.type='button';
    const line=element('div','','line');line.append(element('strong',kindLabel[j.kind]||j.kind),jobStatus(j.status));
    const date=new Date(j.created_at);b.append(line,element('small',String(j.id).slice(0,10)+' · '+(isNaN(date.valueOf())?'':date.toLocaleString('tr-TR'))));
    b.onclick=()=>{state.selected=j.id;state.jobSignature='';showJobs();loadJob(true)};el.append(b);
  }
}
async function loadJobs(force=false){
  if(state.polling)return;state.polling=true;
  try{
    const data=await api('/api/jobs');state.jobs=Array.isArray(data)?data:data.jobs;
    if(state.page==='methods')syncExperimentChoices();
    const count=state.jobs.filter(j=>['running','queued'].includes(j.status)).length;
    $('#queue-count').hidden=!count;$('#queue-count').textContent=count;
    if(state.page==='jobs'){
      if(!state.selected&&state.jobs.length)state.selected=state.jobs[0].id;
      showJobs();await loadJob(force);
    }
  }catch(e){if(state.page==='jobs')notice(e.message,true)}finally{state.polling=false}
}
$('#refresh-jobs').onclick=e=>safe(()=>loadJobs(true),e.currentTarget);
function details(title,value){const box=element('details');box.append(element('summary',title),element('pre',JSON.stringify(value,null,2)));return box}
function action(label,fn){const b=element('button',label,'secondary');b.type='button';b.onclick=()=>safe(fn,b);return b}
async function openStructure(j){
  const data=await api('/api/jobs/'+encodeURIComponent(j.id)+'/structure');
  if(j.kind==='geometry'||j.kind==='packing'){
    f.reset();
    for(const [key,val] of Object.entries(j.config||{})){
      const input=field(key);if(input)input.value=String(val);
    }
    if('acetate_dyad_fraction' in (j.config||{}))field('dyad_measure').value='fraction';
    field('geometry').value=j.config?.geometry||'single';
    field('width_A').value=String(field('geometry').value==='slit'?(j.config?.gap_A??22):(j.config?.water_height_A??22));
    syncBuilder();
  }
  const n=viewer.setStructure(data);state.viewJob=j.id;
  const toggle=$('#toggle-water');
  toggle.hidden=!data.atoms.some(a=>['solvent','water'].includes(a.role));
  toggle.setAttribute('aria-pressed','true');toggle.title='Suyu gizle';$('#water-legend').classList.remove('omitted');
  $('#site-maps').hidden=true;
  if(j.kind==='geometry'||j.kind==='packing'){
    const shape=j.kind==='geometry'?j.result:await api('/api/jobs/'+encodeURIComponent(j.id)+'/files/structure.json');
    $('#site-maps').hidden=!drawFaceMaps($('#face-maps'),data.atoms,shape,formatNum);
  }
  let source=kindLabel[j.kind]||j.kind;
  if(j.kind==='md' && data.label==='last sampled MD frame')source='Son MD karesi · adım '+data.step;
  if(j.kind==='md' && j.config?.profile==='registered')source='Paketlenmiş başlangıç yapısı';
  $('#view-source').textContent=source+' · '+String(j.id).slice(0,10)+' · Å';
  $('#view-detail').textContent=n.toLocaleString('tr-TR')+' atom';go('builder');
  notice('İş dosyasındaki atom koordinatları görüntüleniyor. Ayarları değiştirirseniz tasarım şemasına dönülür.');
}
function summaryLines(result){
  const rows=[];
  if(!result)return rows;
  if(result.assay==='xrd'){
    rows.push(['Ölçüm kaynağı',result.source_doi?'Yayımlanmış ölçüm':'Kullanıcı ölçümü']);
    rows.push(['Islak − kuru polimer etkisi',formatNum(result.summary?.mean)+' nm']);
    for(const r of result.rows||[]){rows.push([r.batch+' · kuru etki',formatNum(r.dry_effect_nm)+' nm'],[r.batch+' · ıslak etki',formatNum(r.wet_effect_nm)+' nm'])}
  }else if(result.assay==='optical'){
    for(const r of result.rows||[])rows.push([r.batch+' · '+r.medium+' · '+r.polymer+' R₂₄',formatNum(r.R24)]);
    for(const [medium,m] of Object.entries(result.summary?.by_medium||{}))rows.push([medium+' · eşlenmiş PCE − PEO',formatNum(m.mean)]);
  }else{
    const labels={n_frames:'Kare sayısı',mean_plane_separation_nm:'Ortalama düzlem aralığı (nm)',contact_cutoff_A:'Temas eşiği (Å)',mean_water_O_contact_fraction:'Temastaki su oksijeni kesri',mean_polymer_heavy_contact_fraction:'Temastaki polimer ağır atomu kesri',sampled_frames:'Örneklenen kare',atom_count:'Atom sayısı',atoms_count:'Atom sayısı',atoms:'Atom sayısı',polymer:'Polimer',water_count:'Su molekülleri',ion_pairs:'NaCl çifti',salt_M_requested:'Hedef NaCl (mol/L)',salt_M_realized:'Gerçekleşen NaCl (mol/L)',nominal_polymer_mass_per_exposed_area_mg_m2:'Geometrik polimer yükü (mg/m²)',seed:'Tohum',accessible_volume_A3:'Geometrik sıvı hacmi (Å³)',center_spacing_A:'Başlangıç tabaka merkezleri arası (Å)',steps:'Adım',duration_ps:'Süre (ps)',final_mobile_temperature_K:'Son mobil grup sıcaklığı (K)',profile:'Profil',engine:'Motor',frames:'Kare',frame_count:'Kare',dp:'Polimerizasyon derecesi',topology_id:'Topoloji kimliği',charge_e:'Toplam yük (e)',chain_count:'Zincir sayısı',minimize_iterations:'Minimizasyon iterasyonu',initial_mobile_fmax_kcal_mol_A:'İlk mobil Fmax (kcal mol⁻¹ Å⁻¹)',post_minimize_mobile_fmax_kcal_mol_A:'Minimizasyon sonrası Fmax (kcal mol⁻¹ Å⁻¹)',final_mobile_fmax_kcal_mol_A:'Son mobil Fmax (kcal mol⁻¹ Å⁻¹)',fitted_mn_g_mol:'Dağılım Mₙ (g/mol)',fitted_dispersity:'Dağılım Đ',sample_mn_g_mol:'Örneklenen Mₙ (g/mol)',sample_dispersity:'Örneklenen Đ'};
    for(const [key,name] of Object.entries(labels))if(result[key]!==null&&['number','string'].includes(typeof result[key]))rows.push([name,formatNum(result[key])]);
  }
  return rows;
}
function showDetail(j){
  const el=$('#job-detail');el.replaceChildren();
  const top=element('div','','status-line');top.append(element('h2',kindLabel[j.kind]||j.kind),jobStatus(j.status));el.append(top,element('p',j.id,'hint'));
  if(j.error)el.append(element('p',typeof j.error==='string'?j.error:JSON.stringify(j.error),'error'));
  const result=j.result;
  if(result){
    for(const [label,value] of summaryLines(result).slice(0,12)){const row=element('div','','tool-row');row.append(element('span',label),element('strong',value));el.append(row)}
  }
  const fs=(j.files||[]).map(x=>typeof x==='string'?x:x.name).filter(Boolean);
  if(result && j.status==='succeeded'){
    const readout=element('div','','readout');el.append(readout);
    if(j.kind==='geometry')drawStructureStats(readout,result,formatNum);
    if(j.kind==='packing' && fs.includes('structure.json')){
      api('/api/jobs/'+encodeURIComponent(j.id)+'/files/structure.json')
        .then(shape=>drawStructureStats(readout,shape,formatNum))
        .catch(error=>readout.append(element('p',error.message,'hint')));
    }
    if(j.kind==='packing_benchmark')drawBenchmark(readout,result,formatNum);
    if(j.kind==='topology' && result.parameter_provenance){
      drawStructureStats(readout,result.parameter_provenance,formatNum);
      readout.append(details('Parametre kaynakları ve atom eşlemesi',{
        MXene_DOI:result.parameter_provenance.native_doi,
        IFF_SHA256:result.parameter_provenance.native_sha256, OPLS_SHA256:result.parameter_provenance.opls_sha256,
        packed_SHA256:result.packed_sha256, topology_SHA256:result.topology_sha256,
        terms:result.term_counts, cutoff_A:result.cutoff_A
      }));
    }
    if(j.kind==='pva_ensemble' && result.members?.length){
      const select=element('select'),label=element('label','Oluşturucuya aktarılacak zincir');
      for(const [i,m] of result.members.entries())select.add(new Option((i+1)+'. zincir · DP'+m.dp+' · '+m.acetate_positions.length+' asetat',i));
      label.append(select);
      readout.append(label,action('Seçili zinciri kullan',()=>{
        const m=result.members[Number(select.value)];
        applyDesign({polymer:'PVA',dp:m.dp,hydrolysis:m.hydrolysis,acetate_pattern:'explicit',acetate_positions:m.acetate_positions,sequence_seed:m.sequence_seed});
        field('acetate_pattern').closest('details').open=true;
        go('builder');notice('Seçilen zincirin DP ve asetat konumları forma aktarıldı.');
      }));
    }
  }
  if(fs.length){
    const buttons=element('p','','row-actions');
    if(fs.some(x=>['packed.json','structure.json'].includes(x)))buttons.append(action('3B yapıyı aç',()=>openStructure(j)));
    if(j.kind==='packing' && j.status==='succeeded')buttons.append(action('Bu yapıdan topoloji üret',async()=>{
      go('methods');await loadRegistered(undefined,j.id);$('#topology-form').scrollIntoView({behavior:'smooth',block:'center'});
    }));
    if(j.kind==='topology' && result?.topology_id)buttons.append(action('Bu topolojiyle dinamik çalıştır',async()=>{
      go('methods');await loadRegistered(j.config.source_job);regTopology.value=result.topology_id;
      if(!regTopology.value)throw Error('Üretilen topoloji kaynak paketlemeyle eşleştirilemedi.');
      reg.elements.minimize_steps.value='25';reg.elements.steps.value='5';
      reg.scrollIntoView({behavior:'smooth',block:'center'});
    }));
    if(j.kind==='packing' && j.status==='succeeded')buttons.append(action('Parametreli MD',async()=>{
      go('methods');await loadRegistered(j.id);$('#registered-form').scrollIntoView({behavior:'smooth',block:'center'});
    }));
    if(j.kind==='md' && j.status==='succeeded' && fs.includes('trajectory-groups.json'))buttons.append(action('Yörüngeyi analiz et',()=>submit('trajectory',{source_job:j.id})));
    if(buttons.childNodes.length)el.append(buttons);
    if(j.kind==='md' && j.status==='succeeded' && fs.some(x=>x.endsWith('.restart') || x==='restart.bin' || x.endsWith('/restart.bin'))){
      const row=element('form','','resume-row'),label=element('label','Ek adım sayısı'),input=element('input'),b=element('button','MD hesabına devam et','secondary');
      input.type='number';input.min='1';input.max='2000';input.value='20';input.required=true;b.type='submit';
      row.onsubmit=e=>{e.preventDefault();safe(()=>submit('md',{restart_job:j.id,md_steps:Number(input.value)}),b)};
      label.append(input);row.append(label,b);el.append(row);
    }
    if(j.kind==='trajectory' && fs.includes('density.csv')){
      const figure=element('figure','','result-plot');el.append(figure);
      densityPlot(j.id,figure).catch(e=>{figure.replaceChildren(element('p',e.message,'hint'))});
    }
    el.append(element('h3','Çıktı dosyaları'));
    const links=element('div','','file-links');for(const f of fs){const a=element('a',f);a.href='/api/jobs/'+encodeURIComponent(j.id)+'/files/'+encodeURIComponent(f);a.download=f;links.append(a)}el.append(links);
  }
  el.append(details('Kayıtlı girdiler',j.config));
  if(result){const concise=Object.fromEntries(Object.entries(result).filter(([k])=>!['atoms','bonds','templates','polymer_bonds'].includes(k)));el.append(details('Sayısal çıktı ve yöntem',concise))}
}
async function densityPlot(id,figure){
  const r=await fetch('/api/jobs/'+encodeURIComponent(id)+'/files/density.csv');
  if(!r.ok)throw Error('Yoğunluk profili okunamadı.');
  const lines=(await r.text()).trim().split(/\r?\n/),columns=lines.shift().split(',');
  const xi=columns.indexOf('mean_normal_position_A'),yi=columns.indexOf('water_O_number_density_A-3');
  if(xi<0||yi<0)throw Error('Yoğunluk profilinde birim ve koordinat alanları eksik.');
  const vals=lines.map(l=>l.split(',')).map(c=>[Number(c[xi]),Number(c[yi])]).filter(p=>p.every(Number.isFinite));
  if(vals.length<2)throw Error('Grafik için iki veya daha fazla konum gerekir.');
  const xmin=Math.min(...vals.map(v=>v[0])),xmax=Math.max(...vals.map(v=>v[0])),ymax=Math.max(...vals.map(v=>v[1]));
  const svg=document.createElementNS(svgNS,'svg');
  svg.setAttribute('viewBox','0 0 510 175');svg.setAttribute('role','img');svg.setAttribute('aria-label','Yüzeye dik su oksijeni sayı yoğunluğu');
  function tag(name,attrs,msg){const x=document.createElementNS(svgNS,name);for(const [k,v] of Object.entries(attrs))x.setAttribute(k,String(v));if(msg)x.textContent=msg;svg.append(x);return x}
  const x=v=>54+425*(v-xmin)/(xmax-xmin||1),y=v=>140-114*v/(ymax||1);
  for(let k=0;k<3;k++){
    const yv=ymax*k/2;tag('line',{x1:54,x2:479,y1:y(yv),y2:y(yv),stroke:'#dbe5df'});
    tag('text',{x:46,y:y(yv)+4,'text-anchor':'end'},formatNum(yv));
  }
  const pts=vals.map(v=>x(v[0])+','+y(v[1])).join(' ');tag('polyline',{points:pts,fill:'none',stroke:'#187769','stroke-width':2});
  tag('text',{x:54,y:161},formatNum(xmin));tag('text',{x:479,y:161,'text-anchor':'end'},formatNum(xmax));tag('text',{x:264,y:165,'text-anchor':'middle'},'Normal konum (Å)');
  tag('text',{x:54,y:14},'Su oksijeni (atom/Å³)');
  figure.replaceChildren(element('figcaption','Yüzeye dik yoğunluk profili'),svg);
}
async function loadJob(force=false){
  if(!state.selected || state.page!=='jobs')return;
  const id=state.selected,j=state.jobs.find(x=>x.id===id);
  const sig=id+'|'+(j?.status||'')+'|'+(j?.updated_at||'');
  if(!force&&sig===state.jobSignature)return;
  const seq=++state.jobRequest;const full=await api('/api/jobs/'+encodeURIComponent(id));
  if(id===state.selected&&seq===state.jobRequest){state.jobSignature=sig;showDetail(full.job||full)}
}
const assays={
  xrd:{labels:{batch:'Parti',bare_dry_nm:'Çıplak kuru (nm)',coated_dry_nm:'Kaplı kuru (nm)',bare_wet_nm:'Çıplak ıslak (nm)',coated_wet_nm:'Kaplı ıslak (nm)'},hint:'Her satıra aynı MXene partisinin dört ölçümünü girin. Sonuç, kaplamanın ıslak ve kuru etkileri arasındaki farktır.'},
  optical:{labels:{batch:'Parti',medium:'Ortam',polymer:'Polimer',A0:'A₀',A24:'A₂₄',B0:'B₀',B24:'B₂₄'},hint:'A örnek absorbansını, B aynı ortam ve polimerin MXene içermeyen arka planını gösterir. Eşlenmiş karşılaştırma için aynı parti ve ortamda PCE ve PEO kodlarını kullanın.'}
};
const assayForm=$('#experiment-form'),assay=()=>assayForm.elements.assay.value;
const assayDraft={xrd:[],optical:[]};
let activeAssay='xrd';
function rowValues(){
  return $$('#assay-grid tbody tr').map(tr=>Object.fromEntries(Array.from(tr.querySelectorAll('input')).map(i=>[i.dataset.key,i.type==='number'?i.valueAsNumber:i.value.trim()])));
}
function addObservation(saved={}){
  const tr=element('tr');
  for(const [key,label] of Object.entries(assays[assay()].labels)){
    const td=element('td'),i=element('input');i.dataset.key=key;i.setAttribute('aria-label',label);i.required=true;
    if(!['batch','medium','polymer'].includes(key)){i.type='number';i.step='any';if(assay()==='xrd')i.min='0.00000001'}else{i.maxLength=120}
    if(saved[key]!==undefined && saved[key]!==null && !Number.isNaN(saved[key]))i.value=String(saved[key]);
    td.append(i);tr.append(td);
  }
  const td=element('td'),del=element('button','×','remove-row');del.type='button';del.setAttribute('aria-label','Satırı kaldır');
  del.onclick=()=>{tr.remove();if(!$('#assay-grid tbody tr'))addObservation()};td.append(del);tr.append(td);$('#assay-grid tbody').append(tr);
}
function initAssay(){
  const tab=element('table'),head=element('thead'),tr=element('tr');for(const label of [...Object.values(assays[assay()].labels),''])tr.append(element('th',label));
  head.append(tr);tab.append(head,element('tbody'));$('#assay-grid').replaceChildren(tab);$('#assay-hint').textContent=assays[assay()].hint;
  if(assayDraft[assay()].length)assayDraft[assay()].forEach(addObservation);else addObservation();
}
assayForm.elements.assay.onchange=()=>{assayDraft[activeAssay]=rowValues();activeAssay=assay();initAssay()};
$('#add-row').onclick=()=>addObservation();
assayForm.onsubmit=e=>{e.preventDefault();safe(()=>submit('experiment',{assay:assay(),rows:rowValues()}),e.submitter)};
$('#reference-run').onclick=e=>safe(()=>submit('experiment',{assay:'xrd',reference:'xing-pei'}),e.currentTarget);
function drawAtlas(){
  const q=$('#record-search').value.toLocaleLowerCase('tr').trim(),m=$('#record-method').value;
  const filtered=state.records.filter(r=>(!m||r.method===m)&&(!q||[r.id,r.polymer,r.condition,r.state,r.note,r.method,r.source_doi,state.sourceNames?.[r.source_doi]].join(' ').toLocaleLowerCase('tr').includes(q)));
  const tab=element('table'),head=element('thead'),hr=element('tr');
  for(const k of ['Polimer','Yöntem','Koşul','Değer','Kaynak'])hr.append(element('th',k));
  head.append(hr);tab.append(head);
  const body=element('tbody');for(const r of filtered){
    const tr=element('tr');tr.append(element('td',r.polymer),element('td',r.method));const cond=element('td',r.condition);
    if(r.note||r.state){const more=element('details');more.append(element('summary','Numune ve not'),element('p',[r.state,r.note].filter(Boolean).join(' · '),'hint'));cond.append(more)}
    tr.append(cond,element('td',formatNum(r.value)+' '+(r.unit||'')));
    const source=element('td');if(r.source_doi){const link=element('a','DOI ↗');link.href='https://doi.org/'+encodeURI(r.source_doi);link.target='_blank';link.rel='noopener noreferrer';link.title=[state.sourceNames?.[r.source_doi],r.source_doi].filter(Boolean).join(' · ');source.append(link)}tr.append(source);body.append(tr);
  }
  tab.append(body);$('#record-table').replaceChildren(tab);$('#record-count').textContent=filtered.length+' / '+state.records.length+' kayıt';
}
async function loadAtlas(){
  if(state.atlasLoaded){drawAtlas();return}
  try{const data=await api('/api/atlas');state.records=data.records||[];state.sourceNames=Object.fromEntries((data.sources||[]).map(s=>[s.doi,s.title]));const selector=$('#record-method');selector.length=1;
    for(const m of [...new Set(state.records.map(r=>r.method))].sort()){selector.add(new Option(m,m))}
    state.atlasLoaded=true;drawAtlas();
  }catch(e){notice(e.message,true)}
}
$('#record-search').oninput=drawAtlas;$('#record-method').onchange=drawAtlas;
const reg=$('#registered-form'), regSource=reg.elements.source_job, regTopology=reg.elements.topology_id;
let registrationRequest=0;
const generator=$('#topology-form'), generatedSource=generator.elements.source_job;
let sourceRequest=0;
async function topologyAvailability() {
  const button=generator.querySelector('button'), help=$('#topology-help');
  const id=generatedSource.value, request=++sourceRequest;
  button.disabled=true;
  if(!id) {
    help.textContent=generatedSource.options.length<=1
      ? 'Tamamlanmış bir paketleme gerekli. Sistem oluştur ekranındaki tuzsuz PVA örneğiyle başlayabilirsiniz.'
      : 'Topoloji modeline uygunluğunu görmek için bir paketleme seçin.';
    return;
  }
  help.textContent='Paketleme bileşimi kontrol ediliyor.';
  try {
    const job=await api('/api/jobs/'+encodeURIComponent(id));
    if(request!==sourceRequest)return;
    const polymer=job.config.polymer||'PEO', pairs=job.result.ion_pairs;
    const reasons=[];
    if(pairs>0)reasons.push(pairs+' NaCl çifti var. Bu su modeli için doğrulanmış NaCl parametreleri bulunmuyor.');
    else if(job.result.salt_M_requested>0)reasons.push('Hedef tuz derişimi sıfırdan büyük. Bu model için 0 M seçerek yeniden paketleyin.');
    if(!['PEO','PVA','bare','PEtOx'].includes(polymer))reasons.push(polymer+' için otomatik topoloji ataması desteklenmiyor.');
    if(polymer==='PEtOx' && (job.config.dp??8)>1)reasons.push('PEtOx birimleri arasındaki torsiyon için kaynak parametre eksik.');
    if(state.iffProfiles?.['iff-water']===false)reasons.push('IFF su/katı kaynağını kurulum kılavuzundaki komutla hazırlayın.');
    button.disabled=reasons.length>0;
    help.textContent=reasons.length?reasons.join(' '):
      polymer+' · '+formatNum(job.result.atom_count)+' atom · '+formatNum(pairs)+' NaCl çifti. İki periyodik kenar da kesimin iki katını aşmalıdır.';
  } catch(error) {
    if(request===sourceRequest)help.textContent='Paketleme bilgisi okunamadı. '+error.message;
  }
}
generatedSource.onchange=topologyAvailability;
generator.onsubmit=e=>{
  e.preventDefault();
  if(generator.reportValidity())safe(()=>submit('topology',{source_job:generatedSource.value,cutoff_A:Number(generator.elements.cutoff_A.value)}),e.submitter);
};
async function loadRegistered(preselect,forBuilder){
  const current=preselect||regSource.value,request=++registrationRequest;
  try{
    const [jr,pr]=await Promise.all([api('/api/jobs'),api('/api/profiles')]);
    if(request!==registrationRequest)return;
    const jobs=(jr.jobs||jr).filter(j=>j.kind==='packing'&&j.status==='succeeded');
    regSource.replaceChildren(new Option('Paketleme seçin',''));
    for(const j of jobs){const label=String(j.id).slice(0,10)+' · '+(j.config?.polymer||new Date(j.created_at).toLocaleString('tr-TR'));regSource.add(new Option(label,j.id))}
    const selected=forBuilder || generatedSource.value;
    generatedSource.replaceChildren(...Array.from(regSource.options,option=>new Option(option.text,option.value)));
    generatedSource.value=jobs.some(j=>j.id===selected)?selected:'';
    topologyAvailability();
    regSource.value=jobs.some(j=>j.id===current)?current:'';
    state.profiles=Array.isArray(pr)?pr:(pr.profiles||[]);await fillTopologies();
  }catch(e){if(request!==registrationRequest)return;reg.querySelector('button').disabled=true;generator.querySelector('button').disabled=true;$('#topology-help').textContent=e.message;$('#reg-help').textContent='Topoloji kayıtları yüklenemedi. '+e.message}
}
let topoRequest=0;
async function fillTopologies(){
  const request=++topoRequest,job=regSource.value;
  regTopology.replaceChildren(new Option('Topoloji seçin',''));reg.querySelector('button').disabled=true;
  if(!job){$('#reg-help').textContent='Paketlenmiş bir sisteme ait kayıtlı parametre takımını seçin.';return}
  try{
    const manifest=await api('/api/jobs/'+encodeURIComponent(job)+'/files/manifest.json');if(request!==topoRequest)return;
    const hash=manifest.files?.['packed.json'];
    for(const p of state.profiles||[]){
      const id=p.topology_id||p.id;
      if(id && hash && p.packed_sha256===hash)regTopology.add(new Option((p.label||p.name||p.profile_name||'Topoloji')+' · '+id.slice(0,12),id));
    }
    if(regTopology.options.length===1){$('#reg-help').textContent='Bu paketlemeye kayıtlı topoloji yok. Yukarıdaki topoloji üreticisini kullanabilir veya kendi parametreli girdinizi Studio CLI ile kaydedebilirsiniz.';return}
    regTopology.selectedIndex=1;reg.querySelector('button').disabled=false;$('#reg-help').textContent='Girdi kimliği eşleşti. Sıfır adım enerji ve kuvvet hesabı üretir.';
  }catch(e){if(request===topoRequest)$('#reg-help').textContent=e.message}
}
regSource.onchange=fillTopologies;
reg.onsubmit=e=>{e.preventDefault();safe(()=>submit('md',{profile:'registered',source_job:regSource.value,topology_id:regTopology.value,steps:Number(reg.elements.steps.value),temperature_k:Number(reg.elements.temperature_k.value),minimize_steps:Number(reg.elements.minimize_steps.value)}),e.submitter)};
const planner=$('#plan-form');
planner.elements.mode.onchange=()=>{$('#ollama-label').hidden=planner.elements.mode.value!=='ollama';planner.elements.model.required=!$('#ollama-label').hidden};
planner.onsubmit=e=>{
  e.preventDefault();const a=planner.elements,body={mode:a.mode.value,objective:a.objective.value};if(body.mode==='ollama')body.model=a.model.value.trim();
  safe(async()=>{
    const v=await api('/api/plan',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const div=$('#plan-result');div.replaceChildren();
    for(const [i,j] of (v.jobs||[]).entries()){
      const card=element('div','','plan-job');card.append(element('strong',kindLabel[j.kind]+' · '+(i+1)),details('Önerilen girdiler',j.config),action('Bu reçeteyi sırala',()=>submit(j.kind,j.config)));div.append(card);
    }
  },e.submitter);
};
const campaign=$('#campaign-form');
try{state.campaign=localStorage.getItem('polychemmd-recent-campaign')}catch{state.campaign=null}
function syncExperimentChoices(){
  const sel=campaign.elements.experiment_id,selected=sel.value,jobs=state.jobs.filter(j=>j.kind==='experiment'&&j.status==='succeeded');
  if(jobs.map(j=>j.id).join(',')===sel.dataset.signature)return;
  sel.dataset.signature=jobs.map(j=>j.id).join(',');sel.replaceChildren(new Option('Bağlantı yok',''));
  for(const j of jobs)sel.add(new Option((j.config?.assay||'Ölçüm')+' · '+String(j.id).slice(0,12),j.id));
  if(jobs.some(j=>j.id===selected))sel.value=selected;
}
campaign.elements.objective.onchange=()=>{$('#campaign-step').hidden=campaign.elements.objective.value!=='model'};
campaign.onsubmit=e=>{
  e.preventDefault();const a=campaign.elements,body={objective:a.objective.value,replicates:Number(a.replicates.value),base_seed:Number(a.base_seed.value)};
  if(body.objective==='model')body.md_steps=Number(a.md_steps.value);if(a.experiment_id.value)body.experiment_id=a.experiment_id.value;
  safe(async()=>{const v=await api('/api/campaigns',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    state.campaign=v.id;try{localStorage.setItem('polychemmd-recent-campaign',v.id)}catch{}
    showCampaign(v);notice('Eşlenmiş seri kaydedildi.');
  },e.submitter);
};
function showCampaign(v){
  const signature=JSON.stringify(v);if(signature===state.campaignSignature)return;state.campaignSignature=signature;
  const el=$('#campaign-result');el.replaceChildren(element('h3','Seri '+String(v.id).slice(0,10)),element('p',(v.complete_pairs??0)+' / '+v.replicates+' tamamlanan çift · '+(statusLabel[v.status]||v.status),'muted'));
  const cols=Array.isArray(v.arms)?v.arms:['Koşul A','Koşul B'],table=element('table','','paired-table'),head=element('thead'),hr=element('tr');
  for(const s of ['Tohum',...cols])hr.append(element('th',typeof s==='string'?s:JSON.stringify(s)));
  head.append(hr);table.append(head);const body=element('tbody');
  for(const p of v.pairs||[]){const row=element('tr');row.append(element('td',String(p.seed)));for(const job of p.jobs||[]){
    const td=element('td'),b=action(statusLabel[job.status]||job.status,()=>{state.selected=job.id;state.jobSignature='';go('jobs');return loadJobs(true)});td.append(b);row.append(td);
  }body.append(row)}
  table.append(body);const wrap=element('div','','table-scroll');wrap.append(table);el.append(wrap);
  if(v.metrics||v.scope)el.append(details('Farklar, birimler ve yöntem',{metrics:v.metrics,units:v.units,scope:v.scope,experiment_id:v.experiment_id}));
}
async function loadCampaign(){
  if(!state.campaign||state.page!=='methods')return;
  try{showCampaign(await api('/api/campaigns/'+encodeURIComponent(state.campaign)))}catch(e){$('#campaign-result').textContent='Son seri bu veri dizininde okunamıyor. '+e.message}
}
setInterval(()=>{if(!document.hidden)loadCampaign()},5500);
initAssay();syncBuilder();go(location.hash.slice(1)||'builder',false);loadStatus();loadJobs();
setInterval(()=>{if(!document.hidden && (state.service || state.page==='jobs'))loadJobs()},3500);
setInterval(()=>{if(!document.hidden)loadStatus()},18000);
