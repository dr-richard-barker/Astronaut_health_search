/* ============================================================================
   Barker Lab shared theme — behaviour (v0.1 prototype)
   - Builds a toggleable left rail with two panels:
       "On this page"  = document map, auto-generated from <section> headings
       "All projects"  = site map, from window.BARKER_SITES (sites.js)
   - Scroll-spy highlights the active section.
   - Remembers open/closed + last panel in localStorage.
   - Progressive enhancement: if this script doesn't run, the page is unchanged.

   To use on a page:  set  <body data-site-id="THIS_REPO_SLUG">
   and include:  theme.css, sites.js, theme.js.  Nothing else required —
   the document map is derived from each <section>'s first <h2>.
   ========================================================================== */
(function(){
  "use strict";
  var LS_OPEN = "barker.map.open", LS_PANEL = "barker.map.panel";
  var slug = document.body.getAttribute("data-site-id") || "";
  var LOGO = document.body.getAttribute("data-brand-logo") || "assets/cose-logo.png";
  var BRAND_URL = document.body.getAttribute("data-brand-url") || "https://cosecloud.com/";

  /* ---------- build DOM shell ---------- */
  var toggle = el("button", {class:"map-toggle", "aria-label":"Toggle navigation map",
    "aria-expanded":"false"});
  toggle.innerHTML = '<span class="bars"><span></span></span><span class="lbl">Map</span>';

  /* top-left bar: COSE brand logo + the toggle */
  var topbar = el("div", {class:"topbar"});
  var brand = el("a", {class:"cose-brand", href:BRAND_URL, target:"_blank",
    rel:"noopener", title:"COSE — cosecloud.com", "aria-label":"COSE — cosecloud.com"});
  brand.appendChild(el("img", {src:LOGO, alt:"COSE"}));
  topbar.appendChild(brand); topbar.appendChild(toggle);

  var bar = el("nav", {class:"sitebar", "aria-label":"Site and document map"});
  var switcher = el("div", {class:"switch", role:"tablist"});
  var tabDoc  = tab("On this page", true);
  var tabSite = tab("All projects", false);
  switcher.appendChild(tabDoc); switcher.appendChild(tabSite);

  var docPanel  = el("div", {class:"panel", id:"panel-doc"});
  var sitePanel = el("div", {class:"panel", id:"panel-site", hidden:""});
  bar.appendChild(switcher); bar.appendChild(docPanel); bar.appendChild(sitePanel);

  var scrim = el("div", {class:"map-scrim"});

  /* wrap existing body content so we can push it when the rail opens */
  var page = el("div", {class:"page"});
  while(document.body.firstChild){ page.appendChild(document.body.firstChild); }
  document.body.appendChild(topbar);
  document.body.appendChild(bar);
  document.body.appendChild(scrim);
  document.body.appendChild(page);

  /* COSE brand link in the footer, immediately before the author name.
     Falls back through common footer shapes (radiation uses footer>.wrap>p;
     the shared template uses a bare <footer class="site">text). */
  var footP = page.querySelector("footer .wrap p") || page.querySelector("footer .wrap")
    || page.querySelector("footer p") || page.querySelector("footer");
  if(footP){
    var fLink = el("a", {class:"cose-foot", href:BRAND_URL, target:"_blank",
      rel:"noopener", title:"COSE — cosecloud.com", "aria-label":"COSE — cosecloud.com"});
    fLink.appendChild(el("img", {src:LOGO, alt:"COSE"}));
    fLink.appendChild(document.createTextNode("COSE"));
    footP.insertBefore(fLink, footP.firstChild);
  }

  /* ---------- document map (auto from <h2> headings) ----------
     Works whether headings are wrapped in <section> (radiation page) or are
     bare <h2 id> children of <main> (the shared Okabe-Ito template). */
  var links = [];
  var scope = page.querySelector("main") || page;
  var heads = [].slice.call(scope.querySelectorAll("h2"));
  var docList = el("ul", {class:"docmap"});
  var idx = 0;
  heads.forEach(function(h){
    if(!h.id){ h.id = "sec-" + (++idx); }
    // Format "<span class=n>01</span>Title" as "01 · Title"; strip inline .tag pills.
    var num = h.querySelector(".n");
    var label;
    if(num){
      var rest = h.textContent.slice(num.textContent.length);
      label = num.textContent.trim() + " · " + rest.replace(/\s+/g," ").trim();
    } else {
      label = h.textContent.replace(/\s+/g," ").trim();
    }
    var pill = h.querySelector(".tag");
    if(pill){ label = label.replace(pill.textContent.replace(/\s+/g," ").trim(),"").trim(); }
    var a = el("a", {href:"#"+h.id});
    a.textContent = label;
    a.addEventListener("click", function(){ if(isOverlay()) close(); });
    var li = el("li"); li.appendChild(a); docList.appendChild(li);
    links.push({a:a, sec:h});
  });
  docPanel.appendChild(header4("On this page"));
  docPanel.appendChild(docList);
  if(!links.length){ docPanel.appendChild(hint("No sections on this page.")); }

  /* ---------- site map (from registry) ---------- */
  sitePanel.appendChild(header4("All projects"));
  var reg = window.BARKER_SITES;
  if(reg && reg.groups){
    var sl = el("ul", {class:"sitemap"});
    reg.groups.forEach(function(g){
      var gl = el("li"); var gh = el("div",{class:"group"}); gh.textContent = g.name;
      gl.appendChild(gh); sl.appendChild(gl);
      g.items.forEach(function(it){
        var li = el("li");
        var a = el("a", {href:it.url});
        if(it.id === slug){ a.className = "current"; }
        a.innerHTML = esc(it.title) + (it.desc? "<small>"+esc(it.desc)+"</small>":"")
          + (it.live===false? '<small style="color:var(--muted)">· page pending</small>':"");
        li.appendChild(a); sl.appendChild(li);
      });
    });
    sitePanel.appendChild(sl);
  } else {
    sitePanel.appendChild(hint("Site registry not loaded."));
  }

  /* ---------- panel switching ---------- */
  function selectPanel(which){
    var doc = which==="doc";
    tabDoc.setAttribute("aria-selected", doc);
    tabSite.setAttribute("aria-selected", !doc);
    docPanel.hidden = !doc; sitePanel.hidden = doc;
    try{ localStorage.setItem(LS_PANEL, which); }catch(e){}
  }
  tabDoc.addEventListener("click", function(){ selectPanel("doc"); });
  tabSite.addEventListener("click", function(){ selectPanel("site"); });

  /* ---------- open / close ---------- */
  function open(){ document.body.classList.add("map-open"); toggle.setAttribute("aria-expanded","true");
    try{ localStorage.setItem(LS_OPEN,"1"); }catch(e){} }
  function close(){ document.body.classList.remove("map-open"); toggle.setAttribute("aria-expanded","false");
    try{ localStorage.setItem(LS_OPEN,"0"); }catch(e){} }
  function isOverlay(){ return window.matchMedia("(max-width:1179px)").matches; }
  toggle.addEventListener("click", function(){
    document.body.classList.contains("map-open") ? close() : open(); });
  scrim.addEventListener("click", close);
  document.addEventListener("keydown", function(e){ if(e.key==="Escape") close(); });

  /* ---------- scroll-spy ---------- */
  var spy = null;
  if("IntersectionObserver" in window && links.length){
    spy = new IntersectionObserver(function(entries){
      entries.forEach(function(en){
        if(en.isIntersecting){
          links.forEach(function(l){ l.a.classList.toggle("active", l.sec===en.target); });
        }
      });
    }, {rootMargin:"-45% 0px -50% 0px", threshold:0});
    links.forEach(function(l){ spy.observe(l.sec); });
  }

  /* ---------- restore state ---------- */
  var savedPanel = "doc";
  try{ savedPanel = localStorage.getItem(LS_PANEL) || "doc"; }catch(e){}
  selectPanel(savedPanel);
  var wantOpen = false;
  try{ wantOpen = localStorage.getItem(LS_OPEN)==="1"; }catch(e){}
  // Default: closed on first visit so the page looks identical to before.
  if(wantOpen) open();

  /* ---------- tiny helpers ---------- */
  function el(tag, attrs){ var n=document.createElement(tag);
    if(attrs) for(var k in attrs){ if(k in n && k!=="hidden" && typeof n[k]!=="object"){} n.setAttribute(k, attrs[k]); }
    return n; }
  function tab(text, sel){ var b=el("button",{role:"tab","aria-selected":String(sel)});
    b.textContent=text; return b; }
  function header4(t){ var h=el("h4"); h.textContent=t; h.style.padding="0 20px"; return h; }
  function hint(t){ var p=el("p"); p.textContent=t;
    p.style.cssText="padding:8px 20px;color:var(--muted);font-size:.85rem"; return p; }
  function esc(s){ return String(s).replace(/[&<>"]/g,function(c){
    return {"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]; }); }
})();
