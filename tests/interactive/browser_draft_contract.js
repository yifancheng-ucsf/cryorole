"use strict";

// Execute the real browser asset with a small deterministic DOM/fetch adapter.
async function runBrowserDraftContracts(source) {
  const results = [];
  async function test(name, action) {
    try { await action(); results.push({name, passed: true}); }
    catch (error) { results.push({name, passed: false, error: String(error.message || error)}); }
  }
  function assert(value, message) { if (!value) throw new Error(message); }
  function immediate(value) { return {then: fn => immediate(fn(value))}; }
  function setup(options = {}) {
    const count = options.count || 3;
    const rv = count === 3 ? [[0,0,0],[0,0,.5],[1,1,.2]] : Array.from({length:count}, (_,i) => [i/count,0,0]);
    const payload = {
      space: "raw", euler_convention: "extrinsic_zyx", run_id: "run", session_token: "session",
      landscape_sha256: "landscape", displayed_point_count: count, full_candidate_count: count,
      sld_field: "sld_display", colormap: "rainbow_r",
      points: {rv, euler:rv, sld:count===3?[1,5,3]:Array.from({length:count},(_,i)=>i+1),
        sld_display_is_outlier:Array(count).fill(false), overlay_selected:Array(count).fill(false)},
    };
    const elements = {};
    const defaults = {representation:"rv", displayMode:"all", filterMode:options.filter||"all", filterValue:options.value||"2", radius:"6", selectionId:"selection"};
    function element(id) {
      if (!elements[id]) {
        const context = {clearRect(){},fillRect(){},beginPath(){},arc(){},fill(){}};
        elements[id] = {id, value:defaults[id]||"", textContent:"", disabled:false, width:620, height:420,
          listeners:{}, addEventListener(type,fn){(this.listeners[type] ||= []).push(fn);},
          getContext(){return context;},getBoundingClientRect(){return {left:0,top:0,width:620,height:420};}};
      }
      return elements[id];
    }
    const document = {getElementById:element,createElement:()=>({click(){}})};
    const pending = [];
    function fetch(path, options) {
      if (path === "/api/session") return immediate({json:()=>payload});
      return new Promise(resolve => pending.push({path, body:JSON.parse(options.body), resolve}));
    }
    const state = new Function("document","fetch","confirm",source+"\nreturn () => ({draft,center});")(document,fetch,()=>true);
    function event(id,type,extra={}) { for(const fn of element(id).listeners[type]||[]) fn({target:element(id),...extra}); }
    function response(index, data, ok=true) { pending[index].resolve({ok,json:async()=>data}); }
    function draft(id, radius=6) {return {draft_id:id,selected_count:1,full_candidate_count:count,selected_fraction:1/count,
      space:"raw",radius,radius_unit:"degrees",center_input_representation:"rotvec",center_input:[0,0,0],
      center_rv:[0,0,0],center_euler:[0,0,0],display_selected:Array(count).fill(false),
      display_local_neighborhood:Array(count).fill(true),command_template:"cryorole select"};}
    async function evaluate(id="current") {
      const index=pending.length, finished=element("evaluate").onclick();response(index,draft(id));await finished;
    }
    return {element,event,response,draft,evaluate,pending,state};
  }
  await test("confirm starts disabled",async()=>{
    const app=setup();assert(app.element("confirm").disabled,"unevaluated Confirm is enabled");
  });
  await test("radius changes invalidate the evaluated draft",async()=>{
    const app=setup();await app.evaluate();app.element("radius").value="7";app.event("radius","input");
    assert(!app.state().draft&&app.element("confirm").disabled,"old radius can still be confirmed");
  });
  await test("point changes invalidate the evaluated draft",async()=>{
    const app=setup();await app.evaluate();app.event("plot0","click",{clientX:600,clientY:20});
    assert(!app.state().draft&&app.element("confirm").disabled,"old center can still be confirmed");
  });
  await test("late evaluation success cannot replace the latest request",async()=>{
    const app=setup(),first=app.element("evaluate").onclick(),second=app.element("evaluate").onclick();
    app.response(1,app.draft("latest"));await second;app.response(0,app.draft("old"));await first;
    assert(app.state().draft.draft_id==="latest","late old result replaced the reviewed draft");
  });
  await test("late evaluation failure cannot replace current feedback",async()=>{
    const app=setup(),first=app.element("evaluate").onclick(),second=app.element("evaluate").onclick();
    app.response(1,app.draft("latest"));await second;app.response(0,{error:"old failure"},false);await first;
    assert(!app.element("message").textContent.includes("old failure"),"stale failure replaced latest success");
  });
  await test("picking ignores hidden particles at the same screen location",async()=>{
    const app=setup({filter:"threshold",value:"2"});app.event("plot0","click",{clientX:35,clientY:390});
    assert(app.state().center[2]===.5,"picked the hidden low-SLD particle");
  });
  await test("top fraction sorts once per redraw",async()=>{
    let count=0;const original=Array.prototype.sort;
    try {
      Array.prototype.sort=function(...args){count++;return original.apply(this,args);};
      setup({count:64,filter:"top",value:"0.4"});
    } finally {Array.prototype.sort=original;}
    assert(count<=1,`sorted the same display ${count} times`);
  });
  await test("large display arrays do not exceed the function argument limit",async()=>{setup({count:140000});});
  await test("confirmation names the exact evaluated snapshot",async()=>{
    const app=setup();await app.evaluate("approved");const index=app.pending.length,finished=app.element("confirm").onclick();
    assert(app.pending[index].body.draft_id==="approved","Confirm did not bind its evaluated snapshot");
    app.response(index,{next_commands:{visualize:"visualize",export:"export"},output_paths:{selection_json:"selection.json"}});await finished;
  });
  await test("display-only changes retain the exact scientific draft",async()=>{
    const app=setup();await app.evaluate("approved");app.element("filterMode").value="threshold";app.event("filterMode","change");
    assert(app.state().draft.draft_id==="approved"&&!app.element("confirm").disabled,"display controls changed the scientific draft");
  });
  return {passed:results.filter(r=>r.passed).length,failed:results.filter(r=>!r.passed).length,results};
}

if(typeof module!=="undefined")module.exports={runBrowserDraftContracts};
if(typeof require!=="undefined"&&require.main===module){
  runBrowserDraftContracts(require("fs").readFileSync(process.argv[2],"utf8")).then(result=>{
    console.log(JSON.stringify(result,null,2));process.exitCode=result.failed?1:0;
  });
}
