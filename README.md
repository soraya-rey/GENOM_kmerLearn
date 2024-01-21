# Classification de gènes par profil de kmer de génomes procaryotes
Ce projet a pour objectif de classifier des gènes de provenance inconnus en utilisant les profil de kmer de génomes.
Nous avons utilisé quatre méthodes:
- Distance entre profils -> contenu dans kmer_distances.ipynb
- Support Vector Machins (SVM) -> contenu dans SVM.ipynb
- Convolutive Neural Network (CCN) entre profiles -> contenu dans classifiers/CNN_with_frequency_matrix/CNN_with_frequency_matrix.ipynb
- CNN entre Chaos Game Representation (CGR) des profiles -> contenu dans CNN_on_CGR.ipynb
Sont également fourni:
- une pipeline de création de la base de données -> database/script.py
- les codes de générations des profiles -> kmer_profiler/generate_profiles.py (et deux autres pour les genes et les protéines)
- les profiles utilisés pour l'entrainement des différentes méthodes