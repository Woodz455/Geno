/** D'où viennent les octets : un tampon en mémoire, ou un serveur qui sait répondre à Range.
 *
 *  Le lecteur ne demande jamais « le fichier » : il demande des plages. C'est toute la
 *  semaine 9 — ne jamais charger ce qu'on ne regarde pas.
 */

export interface Source {
  /** Lit `length` octets à partir de `offset`. Peut rendre moins à la fin du fichier. */
  read(offset: number, length: number): Promise<Uint8Array>;
}

export class BufferSource implements Source {
  constructor(private readonly bytes: Uint8Array) {}

  async read(offset: number, length: number): Promise<Uint8Array> {
    return this.bytes.subarray(offset, Math.min(offset + length, this.bytes.length));
  }
}

export interface RequestRecord {
  offset: number;
  length: number;
  status: number;
  received: number;
  ms: number;
}

/** Requêtes HTTP Range, avec un journal de chacune.
 *
 *  Le journal n'est pas un outil de débogage, c'est l'objet de la mesure : le critère
 *  « premier rendu en une requête » se vérifie en le lisant.
 *
 *  Un serveur qui ignore `Range` répond 200 avec le fichier entier. On l'accepte — on
 *  garde le corps et on sert toutes les lectures suivantes depuis lui — mais on le
 *  note, parce que c'est précisément le comportement que le format cherche à éviter.
 */
export class HttpSource implements Source {
  readonly log: RequestRecord[] = [];
  private whole: Uint8Array | null = null;

  constructor(
    private readonly url: string,
    private readonly fetcher: typeof fetch = (...a) => fetch(...a),
  ) {}

  async read(offset: number, length: number): Promise<Uint8Array> {
    if (this.whole) return this.whole.subarray(offset, offset + length);
    const t0 = performance.now();
    const res = await this.fetcher(this.url, {
      headers: { Range: `bytes=${offset}-${offset + length - 1}` },
      cache: "no-store",
    });
    if (res.status !== 206 && res.status !== 200) {
      throw new Error(`${this.url} : HTTP ${res.status} sur bytes=${offset}-${offset + length - 1}`);
    }
    const body = new Uint8Array(await res.arrayBuffer());
    this.log.push({ offset, length, status: res.status, received: body.length, ms: performance.now() - t0 });
    if (res.status === 200) {
      this.whole = body;
      return body.subarray(offset, offset + length);
    }
    return body;
  }
}
