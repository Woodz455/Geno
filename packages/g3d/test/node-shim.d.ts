// Les deux symboles Node dont le test de lecture croisée a besoin, et rien d'autre.
// `@types/node` entier pour `readFileSync` et `process` serait une dépendance de plus
// pour un fichier de test ; le dépôt n'en a que deux, esbuild et TypeScript.
declare module "node:fs" {
  export function readFileSync(path: string): Uint8Array;
  export function readFileSync(path: string, encoding: "utf8"): string;
}
declare const process: { argv: string[]; exit(code: number): never };
