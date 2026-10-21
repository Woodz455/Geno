// igv.js publie ses types sous `igv`, pas sous le module ES que la page importe (le champ
// `browser` de son package.json désigne la version UMD, qu'esbuild prendrait sinon).
declare module "igv/dist/igv.esm.js" {
  const igv: { createBrowser(div: HTMLElement, options: unknown): Promise<unknown> };
  export default igv;
}
