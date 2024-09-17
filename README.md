# Détection des Répétitions Maximales

## Description
Le projet **Maximal Repeats Detection** est un outil en ligne de commande écrit en Python pour détecter et afficher les répétitions maximales dans une séquence ADN au format FASTA. Une répétition maximale est une séquence qui apparaît au moins deux fois dans l'ADN, mais qui ne peut pas être étendue en raison de la variation des lettres qui l'encadrent.

### Fonctionnalités Principales
- **Détection des Répétitions Maximales** : Identifie toutes les répétitions maximales dans une séquence ADN donnée.
- **Filtrage par Longueur** : Permet de spécifier une longueur minimale et/ou maximale pour les répétitions à rechercher.
- **Exclusion des Répétitions Chevauchantes** : Option pour ignorer les répétitions qui se chevauchent dans la séquence.
- **Export des Résultats** : Les résultats sont enregistrés dans un fichier texte, avec un nom basé sur le script et un numéro incrémenté.

### Définition de la Répétition Maximale
Une **répétition maximale** est une séquence répétée qui ne peut pas être étendue dans la séquence ADN donnée. Pour qu'une répétition soit considérée comme maximale :
- Les occurrences de la séquence doivent être entourées par des lettres différentes autour de chaque occurrence.
- Si une répétition apparaît plus de deux fois, il suffit qu'une paire d'occurrences soit entourée par des lettres différentes pour que la répétition soit considérée comme maximale.

### Exemple de Séquence et Détection
Pour une séquence d'ADN comme celle-ci :
example 
TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTAT

Le script identifiera des répétitions maximales telles que `ATGGGTCCAGAGTTTTGTAATTT` et `TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTAT`, avec les positions et les lettres d'encadrement correspondantes.

## Installation

Assurez-vous que Python 3 est installé sur votre système. Aucun module externe n'est nécessaire pour exécuter ce script.

## Utilisation

1. **Exécution du Script**
   Pour exécuter le script, utilisez la commande suivante dans un terminal :
   ```bash
   python3 maximal_repeats.py <fichier.fasta> [longueur_minimale] [longueur_maximale] [-c]

2. **Exemples de Commandes**

<fichier.fasta> : Le fichier FASTA contenant la séquence ADN à analyser.
[longueur_minimale] : (Optionnel) La longueur minimale des répétitions à rechercher. Si non spécifiée, la longueur minimale par défaut est 20.
[longueur_maximale] : (Optionnel) La longueur maximale des répétitions à rechercher. Si non spécifiée, la recherche inclut toutes les longueurs au-dessus de la longueur minimale.
-c : (Optionnel) Exclut les répétitions chevauchantes.

3. **Arguments**
- Recherche des répétitions de longueur supérieure à 20 :
```bash
   python3 maximal_repeats.py test.fasta
```

- Recherche des répétitions de longueur supérieure à 30 :
```bash
   python3 maximal_repeats.py test.fasta 30
```
- Recherche des répétitions dont la longueur est entre 20 et 30 :
  
```bash
   python3 maximal_repeats.py test.fasta 20 30
```
- Exclure les répétitions chevauchantes avec une longueur minimale de 5 et maximale de 10 :
  
```bash
   python3 maximal_repeats.py test.fasta 5 10 -c
``` 
4. **Résultats** 
Les résultats de l'analyse sont enregistrés dans un fichier texte dans le même répertoire que le script. Les fichiers seront nommés maximal_repeats1.txt, maximal_repeats2.txt, etc. Chaque ligne de ce fichier de sortie contient :

- La séquence répétée
- Sa longueur
- Les positions de chaque occurrence
- Les lettres qui encadrent chaque occurrence
