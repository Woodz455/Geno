# Fixtures — données de test, PAS des données de référence

**Rien ici n'est une source scientifique.** Ces fichiers existent pour une seule raison :
l'environnement de développement n'a pas d'accès réseau aux dépôts de données
(voir [`docs/DATA_SOURCES.md` § 8](../../docs/DATA_SOURCES.md)), et le socle 1D doit
pouvoir être écrit et testé quand même.

## Ce qui est exact

Les **bornes des gènes** sont les coordonnées GRCh38 réelles de ACTB, TP53, BRCA1 et HBB.
Elles servent de repères vérifiables : si le parseur décale d'une base, le test le voit.

## Ce qui est approximatif

Tout le reste — exons, bandes cytogénétiques, cCREs, états chromatiniens, pics CTCF.
Les valeurs sont **plausibles mais fabriquées**. Elles ont la bonne forme, le bon format et
le bon ordre de grandeur, ce qui suffit à exercer la machinerie. Elles n'ont aucune valeur
biologique et ne doivent jamais être citées.

## Garde-fou

Le magasin construit depuis ces fixtures porte `"source": "fixtures"` dans son `index.json`,
et la CLI l'affiche à chaque requête. Une fiche de locus issue des fixtures est reconnaissable
au premier coup d'œil. Quand les vraies données arriveront, `source` portera l'accession et
le hash du manifeste.

## Remplacement

Ces fichiers sont jetables. Dès que `make data` passe, on construit depuis
`pipeline/data/raw/` avec un `tracks.json` pointant sur les vrais fichiers, et les fixtures
ne servent plus qu'aux tests unitaires.
