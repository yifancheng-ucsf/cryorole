"use strict";
let session, draft = null, center = [0, 0, 0], centerRepresentation = "rotvec", sldLow = 0, sldHigh = 1;
const canvases = [0,1,2].map(i => document.getElementById(`plot${i}`));
const pairs = [[0,1],[1,2],[0,2]];

fetch("/api/session").then(r => r.json()).then(data => {
  session = data;
  document.getElementById("spaceBadge").textContent = `${data.space.toUpperCase()} · ${data.euler_convention}`;
  document.getElementById("status").textContent = `Offline localhost · run ${data.run_id}`;
  document.getElementById("summary").textContent = `${data.displayed_point_count.toLocaleString()} displayed / ${data.full_candidate_count.toLocaleString()} full candidates · deterministic display sample · exact selection pending`;
  const outliers=data.points.sld_display_is_outlier.filter(Boolean).length;
  if(data.points.sld.length){sldLow=Math.min(...data.points.sld);sldHigh=Math.max(...data.points.sld);document.getElementById("colorInfo").textContent=`${data.sld_field} · ${data.colormap} · saturated to displayed [${sldLow.toPrecision(4)}, ${sldHigh.toPrecision(4)}] · ${outliers} display outliers`;}else{document.getElementById("colorInfo").textContent=`${data.sld_field} · no points pass the display filter`;}
  draw();
});

function values() { return document.getElementById("representation").value === "euler" ? session.points.euler : session.points.rv; }
function colorForSld(value) {
  const t=Math.max(0,Math.min(1,(value-sldLow)/(sldHigh-sldLow||1)));
  if(session.colormap === "rainbow_r") { const hue=(1-t)*270; return `hsl(${hue} 78% 43%)`; }
  const stops=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]], scaled=t*(stops.length-1), i=Math.min(stops.length-2,Math.floor(scaled)), f=scaled-i;
  return `rgb(${stops[i].map((v,j)=>Math.round(v+(stops[i+1][j]-v)*f)).join(",")})`;
}
function visible(i) {
  const mode = document.getElementById("displayMode").value;
  if (mode === "selected" && draft) return draft.display_selected[i];
  if (mode === "local" && draft) return draft.display_local_neighborhood[i];
  const filter = document.getElementById("filterMode").value, value = Number(document.getElementById("filterValue").value);
  if (filter === "threshold") return session.points.sld[i] >= value;
  if (filter === "top") {
    const sorted = [...session.points.sld].sort((a,b) => b-a), index = Math.max(0, Math.ceil(sorted.length * value) - 1);
    return session.points.sld[i] >= sorted[index];
  }
  return true;
}
function draw() {
  if (!session) return;
  const data = values(), labels = document.getElementById("representation").value === "euler" ? ["α (deg)","β (deg)","γ (deg)"] : ["RV x (rad)","RV y (rad)","RV z (rad)"];
  if(!data.length){canvases.forEach((c,p)=>{c.getContext("2d").clearRect(0,0,c.width,c.height);document.getElementById(`caption${p}`).textContent="No points pass the display filter";});return;}
  pairs.forEach((pair, p) => {
    const c=canvases[p], x=data.map(v=>v[pair[0]]), y=data.map(v=>v[pair[1]]), xmin=Math.min(...x), xmax=Math.max(...x), ymin=Math.min(...y), ymax=Math.max(...y), ctx=c.getContext("2d");
    ctx.clearRect(0,0,c.width,c.height); ctx.fillStyle="#fbfbfc"; ctx.fillRect(0,0,c.width,c.height);
    data.forEach((v,i) => { if(!visible(i)) return; const px=35+(v[pair[0]]-xmin)/(xmax-xmin||1)*(c.width-55), py=c.height-30-(v[pair[1]]-ymin)/(ymax-ymin||1)*(c.height-50); const selected=draft&&draft.display_selected[i]; ctx.fillStyle=selected?"#d00000":(session.points.overlay_selected[i]?"#7b2cbf":colorForSld(session.points.sld[i])); ctx.globalAlpha=session.points.sld_display_is_outlier[i]?0.45:0.78; ctx.beginPath(); ctx.arc(px,py,selected?5:2.5,0,Math.PI*2); ctx.fill(); });
    ctx.globalAlpha=1; document.getElementById(`caption${p}`).textContent=`${labels[pair[0]]} vs ${labels[pair[1]]}`; c._bounds={xmin,xmax,ymin,ymax,pair};
  });
}
canvases.forEach(c => c.addEventListener("click", event => {
  const data=values(); if(!data.length)return; const rect=c.getBoundingClientRect(), bx=(event.clientX-rect.left)/rect.width*c.width, by=(event.clientY-rect.top)/rect.height*c.height, b=c._bounds; let best=0, score=Infinity;
  data.forEach((v,i)=>{ const px=35+(v[b.pair[0]]-b.xmin)/(b.xmax-b.xmin||1)*(c.width-55), py=c.height-30-(v[b.pair[1]]-b.ymin)/(b.ymax-b.ymin||1)*(c.height-50), s=(px-bx)**2+(py-by)**2; if(s<score){score=s;best=i;} });
  center = document.getElementById("representation").value === "euler" ? session.points.euler[best] : session.points.rv[best]; centerRepresentation = document.getElementById("representation").value === "euler" ? "euler" : "rotvec"; document.getElementById("centerCurrent").textContent=`${centerRepresentation}: ${center.map(v=>v.toFixed(5)).join(", ")}`; document.getElementById("message").textContent="Draft center changed. Evaluate to compute the exact full-data result and both representations.";
}));
document.getElementById("radius").addEventListener("input", e => document.getElementById("radiusValue").textContent=`${Number(e.target.value).toFixed(1)}°`);
["representation","displayMode","filterMode","filterValue"].forEach(id => document.getElementById(id).addEventListener("change", draw));
async function post(path, extra) { const response=await fetch(path,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({session_token:session.session_token,run_id:session.run_id,landscape_sha256:session.landscape_sha256,...extra})}); const data=await response.json(); if(!response.ok) throw new Error(data.error||"request failed"); return data; }
document.getElementById("evaluate").onclick=async()=>{ try { draft=await post("/api/evaluate",{center,representation:centerRepresentation,radius_deg:Number(document.getElementById("radius").value)}); document.getElementById("summary").textContent=`EXACT: ${draft.selected_count.toLocaleString()} / ${draft.full_candidate_count.toLocaleString()} selected (${(draft.selected_fraction*100).toFixed(3)}%) · SO(3) geodesic · ${draft.space} · radius ${draft.radius.toFixed(3)} ${draft.radius_unit}`; document.getElementById("centerCurrent").textContent=`${draft.center_input_representation}: ${draft.center_input.map(v=>Number(v).toFixed(5)).join(", ")}`; document.getElementById("centerRv").textContent=draft.center_rv.map(v=>v.toFixed(6)).join(", "); document.getElementById("centerEuler").textContent=draft.center_euler.map(v=>v.toFixed(5)).join(", "); document.getElementById("command").textContent=draft.command_template; document.getElementById("message").textContent="Exact evaluation complete. Review the result before Confirm."; draw(); } catch(e){document.getElementById("message").textContent=e.message;} };
document.getElementById("confirm").onclick=async()=>{ if(!draft){document.getElementById("message").textContent="Evaluate first.";return;} const selection_id=document.getElementById("selectionId").value.trim(); if(!selection_id){document.getElementById("message").textContent="Enter a selection ID; Confirm never overwrites.";return;} if(!confirm(`Create scientific Selection '${selection_id}' from the exact full-data result?`))return; try{const result=await post("/api/confirm",{selection_id}); document.getElementById("command").textContent=`${result.next_commands.visualize}\n${result.next_commands.export}`; document.getElementById("message").textContent=`Created ${result.output_paths.selection_json}. You may visualize/export it, or continue adjusting this draft and Confirm under a new ID.`;}catch(e){document.getElementById("message").textContent=e.message;} };
document.getElementById("downloadDraft").onclick=()=>{if(!draft)return; download("selection_draft.json",JSON.stringify(draft,null,2),"application/json");};
document.getElementById("downloadPreview").onclick=()=>canvases[0].toBlob(blob=>downloadBlob("selection_preview.png",blob));
function download(name,text,type){downloadBlob(name,new Blob([text],{type}));} function downloadBlob(name,blob){const a=document.createElement("a");a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
