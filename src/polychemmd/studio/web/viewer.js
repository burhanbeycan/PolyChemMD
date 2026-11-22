/* Orthographic scientific-coordinate viewer. Chemical bonds are only drawn when
 * present in the input. The initial, scaled illustration is a separate mode. */
const point = (x, y, z) => ({x, y, z});
const palette = {slab:'#75b6a6',pol:'#f3bd7c',wat:'#65afc6',ion:'#c89cde'};
const shade = {Ti:'#74ada2', C:'#a0b4ae', O:'#de9579', F:'#89cabb', H:'#c7d1c9', Na:'#c0a0db', Cl:'#8fa8de'};
function moleculeColor(atom) {
  const role=String(atom.role || '').toLowerCase();
  if (role.includes('poly')) return atom.element==='H' ? '#cbb79a' : palette.pol;
  if (role.includes('water') || role.includes('solvent')) return atom.element==='H' ? '#a4cdd1' : palette.wat;
  return shade[atom.element] || palette.slab;
}
export class StructureView {
  constructor(canvas) {
    this.el=canvas; this.ctx=canvas.getContext('2d'); this.model=null; this.config={geometry:'slit',polymer:'PEO',chains:1};this.pending=false;
    this.reset(); this.schematic(this.config);
    this.resize=new ResizeObserver(()=>this.schedule()); this.resize.observe(canvas);
    canvas.addEventListener('pointerdown',e=>{canvas.setPointerCapture(e.pointerId);this.drag=[e.clientX,e.clientY]});
    canvas.addEventListener('pointermove',e=>{
      if (!this.drag) return;
      this.yaw+=(e.clientX-this.drag[0])*.009;this.pitch=Math.max(-1.47,Math.min(1.47,this.pitch+(e.clientY-this.drag[1])*.009));
      this.drag=[e.clientX,e.clientY];this.schedule();
    });
    for (const name of ['pointerup','pointercancel','lostpointercapture']) canvas.addEventListener(name,()=>{this.drag=null});
    canvas.addEventListener('wheel',e=>{e.preventDefault();this.magnify(Math.exp(-e.deltaY*.0015))},{passive:false});
    canvas.addEventListener('keydown',e=>{
      if (!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','=','-','0'].includes(e.key)) return;
      e.preventDefault();
      if(e.key==='ArrowLeft')this.yaw-=.12;
      if(e.key==='ArrowRight')this.yaw+=.12;
      if(e.key==='ArrowUp')this.pitch=Math.max(-1.47,this.pitch-.12);
      if(e.key==='ArrowDown')this.pitch=Math.min(1.47,this.pitch+.12);
      if(e.key==='+' || e.key==='=')this.magnify(1.12);
      if(e.key==='-')this.magnify(1/1.12);
      if(e.key==='0')this.reset();
      this.schedule();
    });
  }
  reset(){this.yaw=-.46;this.pitch=.40;this.zoom=1;this.schedule()}
  magnify(f){this.zoom=Math.max(.25,Math.min(7,this.zoom*f));this.schedule()}
  schedule(){if(this.pending)return;this.pending=true;requestAnimationFrame(()=>{this.pending=false;this.draw()})}
  schematic(config){this.model=null;this.config=config;this.schedule()}
  setStructure(data){
    if(!Array.isArray(data?.atoms) || !data.atoms.length)throw Error('Koordinat kaydında atom bulunamadı.');
    const atoms=data.atoms;
    if(!atoms.every(a=>[a.x,a.y,a.z].every(Number.isFinite)))throw Error('Yapı geçersiz atom koordinatı içeriyor.');
    const min=[Infinity,Infinity,Infinity],max=[-Infinity,-Infinity,-Infinity];
    for(const a of atoms)for(const [i,k] of ['x','y','z'].entries()){min[i]=Math.min(min[i],a[k]);max[i]=Math.max(max[i],a[k])}
    this.model={atoms,bonds:data.bonds||[],center:min.map((n,i)=>(n+max[i])/2),span:max.map((n,i)=>Math.max(n-min[i],1))};
    this.reset();
    return atoms.length;
  }
  draw(){
    const ctx=this.ctx;if(!ctx)return;
    const {width:w,height:h}=this.el.getBoundingClientRect();if(w<5||h<5)return;
    const dpr=Math.min(2,window.devicePixelRatio||1);
    if(this.el.width!==Math.round(w*dpr)||this.el.height!==Math.round(h*dpr)){this.el.width=Math.round(w*dpr);this.el.height=Math.round(h*dpr)}
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
    const fit=this.model?Math.min(w*(.77/Math.max(1,...this.model.span)),h*(.76/Math.max(1,...this.model.span))):Math.min(w/27,h/16.3);
    const scale=fit*this.zoom,c=Math.cos(this.yaw),s=Math.sin(this.yaw),cp=Math.cos(this.pitch),sp=Math.sin(this.pitch);
    const project=({x,y,z})=>{
      const xx=x*c-y*s, yy=x*s+y*c;
      return [w*.5+xx*scale,h*.49+(yy*sp-z*cp)*scale,yy*cp+z*sp];
    };
    function line(a,b,color,width=1){
      ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();
    }
    if(this.model){
      const {atoms,bonds,center}=this.model;
      const poses=atoms.map(a=>project(point(a.x-center[0],a.y-center[1],a.z-center[2])));
      for(const b of bonds){
        if(!Array.isArray(b) || !poses[b[0]] || !poses[b[1]])continue;
        line(poses[b[0]],poses[b[1]],'#87988e',.8);
      }
      const order=Array.from({length:atoms.length},(_,i)=>i).sort((a,b)=>poses[a][2]-poses[b][2]);
      for(const i of order){
        const a=atoms[i],p=poses[i],factor=a.element==='H'?.18:a.element==='Ti'?.43:.32;
        const rad=Math.max(.8,Math.min(9,scale*factor));
        ctx.fillStyle=moleculeColor(a);ctx.beginPath();ctx.arc(p[0],p[1],rad,0,2*Math.PI);ctx.fill();
      }
    }else{
      const dual=this.config.geometry==='slit';
      const floor=-3.15,roof=3.15;
      for(let i=-5;i<=5;i++){line(project(point(i*2.2,-7,floor-.55)),project(point(i*2.2,7,floor-.55)),'#29444a',.55);line(project(point(-11,i*1.35,floor-.55)),project(point(11,i*1.35,floor-.55)),'#29444a',.55)}
      const parts=[];
      function layer(z) {
        const ax=-8.85,ay=-3.45;
        const pts=[];for(let j=0;j<=7;j++){const row=[];for(let i=0;i<=17;i++)row.push(project(point(ax+i*1.04+(j%2)*.36,ay+j*.98,z)));pts.push(row)}
        for(let j=0;j<7;j++)for(let i=0;i<17;i++){
          const poly=[pts[j][i],pts[j][i+1],pts[j+1][i+1],pts[j+1][i]];
          parts.push({z:poly.reduce((a,p)=>a+p[2],0)/4,draw:()=>{
            ctx.beginPath();ctx.moveTo(...poly[0].slice(0,2));for(const p of poly.slice(1))ctx.lineTo(p[0],p[1]);ctx.closePath();
            ctx.fillStyle='#285855';ctx.fill();ctx.strokeStyle='#72a89a';ctx.lineWidth=.48;ctx.stroke();line(poly[0],poly[2],'#67988c',.43);
            ctx.fillStyle='#a7c9b8';ctx.beginPath();ctx.arc(poly[0][0],poly[0][1],1.06,0,Math.PI*2);ctx.fill();
          }});
        }
      }
      layer(floor);if(dual)layer(roof);
      const withPoly=this.config.polymer!=='bare' && Number(this.config.chains)!==0;
      if(withPoly){
        let old=null;for(let i=0;i<22;i++){const a=project(point(-7+i*.67,-.85+.38*(i%2)+Math.sin(i*.6)*.48,-.3+Math.sin(i*.57)*.72+.045*i));
          if(old){const b=old;parts.push({z:(a[2]+b[2])/2,draw:()=>{line(a,b,'#e7b477',Math.max(1.4,w/205));ctx.fillStyle='#f7cc95';ctx.beginPath();ctx.arc(a[0],a[1],Math.max(1.6,w/190),0,2*Math.PI);ctx.fill()}})}old=a;
        }
      }
      for(let i=0;i<28;i++){
        const p=project(point(Math.sin(i*14.51)*8,Math.cos(i*31.73)*2.8,Math.sin(i*3.52)*1.65));
        parts.push({z:p[2],draw:()=>{ctx.beginPath();ctx.fillStyle='#629eaf';ctx.arc(p[0],p[1],Math.max(1.4,w/265),0,Math.PI*2);ctx.fill()}});
      }
      parts.sort((a,b)=>a.z-b.z).forEach(x=>x.draw());
    }
    const origin=[25,h-31,0],center=project(point(0,0,0));
    for(const [label,vec] of [['x',point(1,0,0)],['y',point(0,1,0)],['z',point(0,0,1)]]){
      const v=project(vec),n=[origin[0]+(v[0]-center[0])*15/scale,origin[1]+(v[1]-center[1])*15/scale];
      line(origin,n,'#9dbbb3',1);ctx.fillStyle='#aec7bf';ctx.font='9px system-ui';ctx.fillText(label,n[0]+1,n[1]-2);
    }
  }
}
