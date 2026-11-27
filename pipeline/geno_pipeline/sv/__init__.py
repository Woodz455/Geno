"""Mode variant structural (semaine 15) : supprimer une frontière, voir l'enhancer changer de cible.

Le cas est celui de Lupiáñez 2015 (Cell 161:1012) : au locus WNT6/IHH/EPHA4/PAX3, des
délétions de 1,75 à 1,9 Mb emportent EPHA4 et la frontière entre son TAD et celui de PAX3 ;
les enhancers de membre du TAD d'EPHA4 touchent alors PAX3, qui s'exprime à tort dans la
main — une brachydactylie. Une délétion de même forme qui **laisse** la frontière (« DelBS »
chez la souris) ne fait rien.

- `locus` — le locus en hg38 : gènes GENCODE, sites CTCF d'ENCODE (GM12878, orientés par
  motif), enhancers et variants **reconstruits depuis le texte de l'article** (sa Table S1 est
  inaccessible), avec la plage de points de cassure que le texte autorise ;
- `predict` — ce que l'extrusion 1D seule dit des contacts : une distance génomique effective,
  raccourcie par les boucles, convertie en probabilité de contact par une loi calée sur le 3D
  du type sauvage, puis éprouvée sur le 3D des allèles mutants.
"""
