/* Semaine 22 — compteurs WebGL2, injectés avant la page (Page.addScriptToEvaluateOnNewDocument).
 *
 * Ils comptent ce que le moteur demande au GPU, sans toucher au moteur ni dépendre de la
 * machine : appels de dessin, instances, sommets, pixels des passes plein écran, octets alloués.
 * Une image se termine au premier dessin dans le tampon par défaut (la composition finale).
 * Résultat : window.__glcount = { frames: [{ draws, instances, vertices, fsPixels }], bytes, peak }.
 */
(() => {
  const P = WebGL2RenderingContext.prototype;
  const G = WebGL2RenderingContext;
  const BPP = new Map([
    [G.RGBA8, 4], [G.RGB8, 3], [G.RG8, 2], [G.R8, 1], [G.SRGB8_ALPHA8, 4], [G.RGBA, 4], [G.RGB, 3],
    [G.LUMINANCE, 1], [G.ALPHA, 1], [G.LUMINANCE_ALPHA, 2],
    [G.RGBA16F, 8], [G.RG16F, 4], [G.R16F, 2], [G.RGBA32F, 16], [G.RG32F, 8], [G.R32F, 4], [G.R11F_G11F_B10F, 4],
    [G.RGBA32UI, 16], [G.RG32UI, 8], [G.R32UI, 4], [G.RGBA16UI, 8], [G.RGBA8UI, 4], [G.R8UI, 1], [G.R16UI, 2],
    [G.RGBA32I, 16], [G.R32I, 4],
    [G.DEPTH_COMPONENT16, 2], [G.DEPTH_COMPONENT24, 4], [G.DEPTH_COMPONENT32F, 4],
    [G.DEPTH24_STENCIL8, 4], [G.DEPTH32F_STENCIL8, 8], [G.DEPTH_COMPONENT, 4], [G.DEPTH_STENCIL, 4],
  ]);
  const BIND = new Map([
    [G.ARRAY_BUFFER, G.ARRAY_BUFFER_BINDING], [G.ELEMENT_ARRAY_BUFFER, G.ELEMENT_ARRAY_BUFFER_BINDING],
    [G.UNIFORM_BUFFER, G.UNIFORM_BUFFER_BINDING], [G.COPY_READ_BUFFER, G.COPY_READ_BUFFER_BINDING],
    [G.COPY_WRITE_BUFFER, G.COPY_WRITE_BUFFER_BINDING], [G.PIXEL_PACK_BUFFER, G.PIXEL_PACK_BUFFER_BINDING],
    [G.PIXEL_UNPACK_BUFFER, G.PIXEL_UNPACK_BUFFER_BINDING], [G.TRANSFORM_FEEDBACK_BUFFER, G.TRANSFORM_FEEDBACK_BUFFER_BINDING],
  ]);
  const TEXBIND = new Map([[G.TEXTURE_2D, G.TEXTURE_BINDING_2D], [G.TEXTURE_2D_ARRAY, G.TEXTURE_BINDING_2D_ARRAY], [G.TEXTURE_3D, G.TEXTURE_BINDING_3D]]);
  const sizes = new Map();          // objet GL → octets
  const levels = new Map();         // texture → (niveau → octets), pour texImage2D niveau par niveau
  const state = { frames: [], cur: null, bytes: 0, peak: 0, unknownFormats: 0, fb: null, vp: [0, 0, 0, 0] };
  window.__glcount = state;
  const total = () => {
    let t = 0;
    for (const v of sizes.values()) t += v;
    state.bytes = t;
    state.peak = Math.max(state.peak, t);
  };
  const bpp = (fmt) => {
    const b = BPP.get(fmt);
    if (b === undefined) state.unknownFormats++;
    return b ?? 4;
  };
  const frame = () => (state.cur ??= { draws: 0, instances: 0, vertices: 0, fsPixels: 0, ...(window.__glcountDetail ? { list: [] } : {}) });
  const note = (f, what) => { if (f.list) f.list.push(what); };
  const endIfDefault = () => {
    if (state.fb === null) {
      state.frames.push(state.cur);
      state.cur = null;
    }
  };
  const wrap = (name, fn) => {
    const orig = P[name];
    if (!orig) return;
    P[name] = function (...a) {
      const r = orig.apply(this, a);
      fn.call(this, a);
      return r;
    };
  };

  wrap("bindFramebuffer", function (a) {
    if (a[0] === this.FRAMEBUFFER || a[0] === this.DRAW_FRAMEBUFFER) state.fb = a[1];
  });
  wrap("viewport", (a) => { state.vp = a; });
  wrap("drawArrays", function (a) {
    const f = frame();
    f.draws++; f.instances++; f.vertices += a[2];
    if (a[0] === this.TRIANGLES && a[2] === 3) f.fsPixels += state.vp[2] * state.vp[3];
    note(f, ["arrays", a[2], 1, state.fb ? 1 : 0]);
    endIfDefault();
  });
  wrap("drawElements", function (a) {
    const f = frame();
    f.draws++; f.instances++; f.vertices += a[1];
    endIfDefault();
  });
  wrap("drawArraysInstanced", function (a) {
    const f = frame();
    f.draws++; f.instances += a[3]; f.vertices += a[2] * a[3];
    note(f, ["arraysInstanced", a[2], a[3], state.fb ? 1 : 0]);
    endIfDefault();
  });
  wrap("drawElementsInstanced", function (a) {
    const f = frame();
    f.draws++; f.instances += a[4]; f.vertices += a[1] * a[4];
    note(f, ["elementsInstanced", a[1], a[4], state.fb ? 1 : 0]);
    endIfDefault();
  });

  wrap("bufferData", function (a) {
    const buf = this.getParameter(BIND.get(a[0]));
    if (!buf) return;
    const src = a[1];
    let n = typeof src === "number" ? src : src ? src.byteLength : 0;
    if (typeof src !== "number" && src && a.length >= 4) {     // (target, src, usage, srcOffset[, length])
      const el = src.BYTES_PER_ELEMENT ?? 1;
      n = a.length >= 5 && a[4] ? a[4] * el : n - a[3] * el;
    }
    sizes.set(buf, n);
    total();
  });
  wrap("texStorage2D", function (a) {
    const tex = this.getParameter(TEXBIND.get(a[0]) ?? this.TEXTURE_BINDING_2D);
    let n = 0, w = a[3], h = a[4];
    for (let l = 0; l < a[1]; l++) { n += w * h * bpp(a[2]); w = Math.max(1, w >> 1); h = Math.max(1, h >> 1); }
    sizes.set(tex, n);
    total();
  });
  wrap("texStorage3D", function (a) {
    const tex = this.getParameter(TEXBIND.get(a[0]) ?? this.TEXTURE_BINDING_3D);
    let n = 0, w = a[3], h = a[4];
    for (let l = 0; l < a[1]; l++) { n += w * h * a[5] * bpp(a[2]); w = Math.max(1, w >> 1); h = Math.max(1, h >> 1); }
    sizes.set(tex, n);
    total();
  });
  wrap("texImage2D", function (a) {
    const tex = this.getParameter(TEXBIND.get(a[0]) ?? this.TEXTURE_BINDING_2D);
    if (!tex) return;
    let w, h;
    if (a.length >= 8) { w = a[3]; h = a[4]; } else { w = a[5]?.width ?? 0; h = a[5]?.height ?? 0; }
    const per = levels.get(tex) ?? new Map();
    per.set(a[1], w * h * bpp(a[2]));
    levels.set(tex, per);
    let n = 0;
    for (const v of per.values()) n += v;
    sizes.set(tex, n);
    total();
  });
  wrap("texImage3D", function (a) {
    const tex = this.getParameter(TEXBIND.get(a[0]) ?? this.TEXTURE_BINDING_3D);
    if (!tex) return;
    sizes.set(tex, a[3] * a[4] * a[5] * bpp(a[2]));
    total();
  });
  wrap("renderbufferStorage", function (a) {
    sizes.set(this.getParameter(this.RENDERBUFFER_BINDING), a[2] * a[3] * bpp(a[1]));
    total();
  });
  wrap("renderbufferStorageMultisample", function (a) {
    sizes.set(this.getParameter(this.RENDERBUFFER_BINDING), a[1] * a[3] * a[4] * bpp(a[2]));
    total();
  });
  for (const name of ["deleteBuffer", "deleteTexture", "deleteRenderbuffer"]) {
    wrap(name, (a) => { sizes.delete(a[0]); levels.delete(a[0]); total(); });
  }
})();
