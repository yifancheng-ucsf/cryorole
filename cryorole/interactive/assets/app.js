"use strict";
let session, draft = null, draftRevision = 0, center = [0, 0, 0], centerRepresentation = "rotvec", sldLow = 0, sldHigh = 1;
const canvases = [0,1,2].map(i => document.getElementById(`plot${i}`));
const pairs = [[0,1],[1,2],[0,2]];
const confirmButton = document.getElementById("confirm");
confirmButton.disabled = true;
document.getElementById("evaluate").disabled = true;

function valueRange(values) {
  let low = Infinity, high = -Infinity;
  for (const value of values) { low = Math.min(low, value); high = Math.max(high, value); }
  return [low, high];
}

fetch("/api/session").then(r => r.json()).then(data => {
  session = data;
  document.getElementById("evaluate").disabled = false;
  document.getElementById("spaceBadge").textContent = `${data.space.toUpperCase()} · ${data.euler_convention}`;
  document.getElementById("status").textContent = `Offline localhost · run ${data.run_id}`;
  document.getElementById("summary").textContent = `${data.displayed_point_count.toLocaleString()} displayed / ${data.full_candidate_count.toLocaleString()} full candidates · deterministic display sample · exact selection pending`;
  const outliers = data.points.sld_display_is_outlier.filter(Boolean).length;
  if (data.points.sld.length) {
    [sldLow, sldHigh] = valueRange(data.points.sld);
    document.getElementById("colorInfo").textContent = `${data.sld_field} · ${data.colormap} · saturated to displayed [${sldLow.toPrecision(4)}, ${sldHigh.toPrecision(4)}] · ${outliers} display outliers`;
  } else document.getElementById("colorInfo").textContent = `${data.sld_field} · no points pass the display filter`;
  draw();
});

function values() { return !session ? [] : (document.getElementById("representation").value === "euler" ? session.points.euler : session.points.rv); }
function colorForSld(value) {
  const t = Math.max(0, Math.min(1, (value-sldLow)/(sldHigh-sldLow||1)));
  if (session.colormap === "rainbow_r") { const hue=(1-t)*270; return `hsl(${hue} 78% 43%)`; }
  const stops=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]], scaled=t*(stops.length-1), i=Math.min(stops.length-2,Math.floor(scaled)), f=scaled-i;
  return `rgb(${stops[i].map((v,j)=>Math.round(v+(stops[i+1][j]-v)*f)).join(",")})`;
}

function visibleMask() {
  const mode = document.getElementById("displayMode").value;
  if (draft && mode === "selected") return draft.display_selected;
  if (draft && mode === "local") return draft.display_local_neighborhood;
  const filter = document.getElementById("filterMode").value, value = Number(document.getElementById("filterValue").value);
  if (filter === "threshold") return session.points.sld.map(sld => sld >= value);
  if (filter === "top") {
    if (!(value > 0 && value <= 1)) return session.points.sld.map(() => false);
    const sorted = [...session.points.sld].sort((a,b) => b-a);
    const cutoff = sorted[Math.max(0, Math.ceil(sorted.length * value)-1)];
    return session.points.sld.map(sld => sld >= cutoff);
  }
  return session.points.sld.map(() => true);
}

function draw() {
  if (!session) return;
  const data = values(), mask = visibleMask();
  const labels = document.getElementById("representation").value === "euler" ? ["α (deg)","β (deg)","γ (deg)"] : ["RV x (rad)","RV y (rad)","RV z (rad)"];
  if (!data.length) {
    canvases.forEach((c,p) => { c.getContext("2d").clearRect(0,0,c.width,c.height); document.getElementById(`caption${p}`).textContent="No points pass the display filter"; });
    return;
  }
  const bounds = [[Infinity,-Infinity],[Infinity,-Infinity],[Infinity,-Infinity]];
  for (const row of data) for (let axis=0; axis<3; axis++) {
    bounds[axis][0] = Math.min(bounds[axis][0],row[axis]);
    bounds[axis][1] = Math.max(bounds[axis][1],row[axis]);
  }
  pairs.forEach((pair,p) => {
    const c=canvases[p], [xmin,xmax]=bounds[pair[0]], [ymin,ymax]=bounds[pair[1]], ctx=c.getContext("2d");
    ctx.clearRect(0,0,c.width,c.height); ctx.fillStyle="#fbfbfc"; ctx.fillRect(0,0,c.width,c.height);
    data.forEach((v,i) => {
      if (!mask[i]) return;
      const px=35+(v[pair[0]]-xmin)/(xmax-xmin||1)*(c.width-55), py=c.height-30-(v[pair[1]]-ymin)/(ymax-ymin||1)*(c.height-50);
      const selected=draft&&draft.display_selected[i];
      ctx.fillStyle=selected?"#d00000":(session.points.overlay_selected[i]?"#7b2cbf":colorForSld(session.points.sld[i]));
      ctx.globalAlpha=session.points.sld_display_is_outlier[i]?0.45:0.78;
      ctx.beginPath(); ctx.arc(px,py,selected?5:2.5,0,Math.PI*2); ctx.fill();
    });
    ctx.globalAlpha=1; document.getElementById(`caption${p}`).textContent=`${labels[pair[0]]} vs ${labels[pair[1]]}`;
    c._bounds={xmin,xmax,ymin,ymax,pair};
  });
}

function invalidateDraft(message) {
  draftRevision++;
  draft=null;
  confirmButton.disabled=true;
  document.getElementById("summary").textContent="Exact evaluation needed for the current center and radius.";
  document.getElementById("centerRv").textContent="—";
  document.getElementById("centerEuler").textContent="—";
  document.getElementById("command").textContent="Evaluate the current center and radius before Confirm.";
  document.getElementById("message").textContent=message;
}

canvases.forEach(c => c.addEventListener("click", event => {
  const data=values(); if (!data.length || !c._bounds) return;
  const mask=visibleMask(), rect=c.getBoundingClientRect(), bx=(event.clientX-rect.left)/rect.width*c.width, by=(event.clientY-rect.top)/rect.height*c.height, b=c._bounds;
  let best=-1, score=Infinity;
  data.forEach((v,i) => {
    if (!mask[i]) return;
    const px=35+(v[b.pair[0]]-b.xmin)/(b.xmax-b.xmin||1)*(c.width-55), py=c.height-30-(v[b.pair[1]]-b.ymin)/(b.ymax-b.ymin||1)*(c.height-50), distance=(px-bx)**2+(py-by)**2;
    if (distance<score) { score=distance; best=i; }
  });
  if (best<0) return;
  centerRepresentation=document.getElementById("representation").value === "euler" ? "euler" : "rotvec";
  center=(centerRepresentation === "euler" ? session.points.euler[best] : session.points.rv[best]).slice();
  invalidateDraft("Draft center changed. Evaluate to compute the exact full-data result and both representations.");
  document.getElementById("centerCurrent").textContent=`${centerRepresentation}: ${center.map(v=>v.toFixed(5)).join(", ")}`;
  draw();
}));

document.getElementById("radius").addEventListener("input", event => {
  document.getElementById("radiusValue").textContent=`${Number(event.target.value).toFixed(1)}°`;
  invalidateDraft("Radius changed. Evaluate the current draft before Confirm.");
  draw();
});
["representation","displayMode","filterMode","filterValue"].forEach(id => document.getElementById(id).addEventListener("change", () => {
  if (id === "filterMode" && document.getElementById(id).value === "top") {
    const input=document.getElementById("filterValue"), value=Number(input.value);
    if (!(value>0 && value<=1)) input.value="0.4";
  }
  draw();
}));

async function post(path, extra) {
  const response=await fetch(path,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({session_token:session.session_token,run_id:session.run_id,landscape_sha256:session.landscape_sha256,...extra})});
  const data=await response.json(); if (!response.ok) throw new Error(data.error||"request failed"); return data;
}

document.getElementById("evaluate").onclick=async()=>{
  invalidateDraft("Evaluating the exact full-data selection…");
  const revision=draftRevision;
  const request={center:center.slice(),representation:centerRepresentation,radius_deg:Number(document.getElementById("radius").value)};
  try {
    const result=await post("/api/evaluate",request);
    if (revision!==draftRevision) return;
    draft=result; confirmButton.disabled=false;
    document.getElementById("summary").textContent=`EXACT: ${draft.selected_count.toLocaleString()} / ${draft.full_candidate_count.toLocaleString()} selected (${(draft.selected_fraction*100).toFixed(3)}%) · SO(3) geodesic · ${draft.space} · radius ${draft.radius.toFixed(3)} ${draft.radius_unit}`;
    document.getElementById("centerCurrent").textContent=`${draft.center_input_representation}: ${draft.center_input.map(v=>Number(v).toFixed(5)).join(", ")}`;
    document.getElementById("centerRv").textContent=draft.center_rv.map(v=>v.toFixed(6)).join(", ");
    document.getElementById("centerEuler").textContent=draft.center_euler.map(v=>v.toFixed(5)).join(", ");
    document.getElementById("command").textContent=draft.command_template;
    document.getElementById("message").textContent="Exact evaluation complete. Review the result before Confirm.";
    draw();
  } catch(error) { if (revision===draftRevision) document.getElementById("message").textContent=error.message; }
};

confirmButton.onclick=async()=>{
  const approved=draft;
  if (!approved) { document.getElementById("message").textContent="Evaluate first."; return; }
  const selection_id=document.getElementById("selectionId").value.trim();
  if (!selection_id) { document.getElementById("message").textContent="Enter a selection ID; Confirm never overwrites."; return; }
  if (!confirm(`Create scientific Selection '${selection_id}' from the exact full-data result (${approved.radius} ${approved.radius_unit})?`)) return;
  confirmButton.disabled=true;
  try {
    const result=await post("/api/confirm",{selection_id,draft_id:approved.draft_id});
    if (draft===approved) document.getElementById("command").textContent=`${result.next_commands.visualize}\n${result.next_commands.export}`;
    document.getElementById("message").textContent=`Created ${result.output_paths.selection_json}. Evaluate any changed center or radius before another Confirm.`;
  } catch(error) { document.getElementById("message").textContent=error.message; }
  finally { confirmButton.disabled=!draft; }
};
document.getElementById("downloadDraft").onclick=()=>{if(draft)download("selection_draft.json",JSON.stringify(draft,null,2),"application/json");};
document.getElementById("downloadPreview").onclick=()=>canvases[0].toBlob(blob=>downloadBlob("selection_preview.png",blob));
function download(name,text,type){downloadBlob(name,new Blob([text],{type}));}
function downloadBlob(name,blob){const a=document.createElement("a");a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
