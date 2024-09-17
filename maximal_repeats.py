import argparse  # Module pour gérer les arguments en ligne de commande
import os  # Module pour gérer les fichiers et les chemins
from Bio import SeqIO  # Module pour lire des fichiers FASTA, une partie de Biopython
       

# Fonction pour lire un fichier FASTA et retourner l'en-tête et la séquence ADN
def parse_fasta(fasta_file):
    """Lit le fichier FASTA et retourne l'en-tête et la séquence."""
    for record in SeqIO.parse(fasta_file, "fasta"):  # Parcourt chaque séquence dans le fichier FASTA
        return record.id, str(record.seq)  # Retourne l'identifiant de la séquence et sa séquence ADN

# Fonction principale pour trouver les répétitions maximales dans la séquence ADN
def find_repeats(sequence, min_len=20, max_len=None, no_overlap=False):
    """Trouve et retourne les répétitions maximales dans une séquence ADN."""
    repeats = []  # Liste pour stocker les répétitions trouvées
    seq_len = len(sequence)  # Longueur de la séquence ADN
    # Si aucune longueur maximale n'est donnée, utiliser la longueur totale de la séquence
    max_len = max_len or seq_len

    # Parcours de toutes les sous-séquences possibles dans la plage de longueurs définie
    for length in range(min_len, max_len + 1):  # Parcourt toutes les longueurs de sous-séquences
        for start in range(seq_len - length + 1):  # Pour chaque position dans la séquence
            word = sequence[start:start + length]  # Extrait le "mot" de longueur donnée
            # Trouve toutes les positions où ce mot apparaît dans la séquence
            positions = [i for i in range(seq_len - length + 1) if sequence[i:i + length] == word]

            # Si le mot apparaît plus d'une fois, il pourrait être une répétition
            if len(positions) > 1:
                encodings = get_encodings(sequence, positions, length)  # Récupère les lettres qui encadrent chaque répétition
                if is_maximal_repeat(encodings):  # Vérifie si c'est une répétition maximale
                    repeats.append((word, length, positions, encodings))  # Ajoute la répétition maximale à la liste

    # Si l'option pour exclure les répétitions chevauchantes est activée, les filtrer
    if no_overlap:
        repeats = filter_non_overlapping(repeats)

    return repeats  # Retourne la liste des répétitions maximales

# Fonction pour obtenir les lettres encadrant chaque occurrence d'une répétition
def get_encodings(sequence, positions, length):
    """Retourne les lettres qui encadrent chaque occurrence d'une répétition."""
    encodings = []  # Liste pour stocker les lettres qui encadrent chaque répétition
    for pos in positions:  # Parcourt chaque position d'occurrence de la répétition
        left = sequence[pos - 1] if pos > 0 else None  # Lettre à gauche (None si en début de séquence)
        right = sequence[pos + length] if pos + length < len(sequence) else None  # Lettre à droite (None si en fin de séquence)
        encodings.append((left, right))  # Stocke le couple (lettre gauche, lettre droite)
    return encodings  # Retourne la liste des encadrements

# Fonction pour vérifier si une répétition est maximale
def is_maximal_repeat(encodings):
    """Vérifie si une répétition est maximale."""
    # Vérifie si les lettres encadrantes sont différentes entre au moins deux occurrences
    return any(encodings[i] != encodings[j] for i in range(len(encodings)) for j in range(i + 1, len(encodings)))

# Fonction pour filtrer les répétitions chevauchantes
def filter_non_overlapping(repeats):
    """Filtre les répétitions chevauchantes."""
    non_overlapping = []  # Liste pour stocker les répétitions non chevauchantes
    for repeat in repeats:  # Parcourt chaque répétition trouvée
        positions = repeat[2]  # Récupère les positions de la répétition
        # Vérifie que les positions ne se chevauchent pas
        if all(positions[i] + repeat[1] <= positions[i + 1] for i in range(len(positions) - 1)):
            non_overlapping.append(repeat)  # Ajoute les répétitions non chevauchantes à la liste
    return non_overlapping  # Retourne la liste des répétitions non chevauchantes

# Fonction pour générer un nom de fichier de sortie qui s'incrémente automatiquement
def generate_output_filename(script_name):
    """Génère un nom de fichier de sortie en fonction du nom du script."""
    base_name = os.path.splitext(script_name)[0]  # Retire l'extension .py du nom du script
    counter = 1  # Initialisation d'un compteur pour générer un numéro de fichier
    output_filename = f"{base_name}{counter}.txt"  # Premier nom de fichier généré

    # Boucle jusqu'à trouver un nom de fichier qui n'existe pas encore
    while os.path.exists(output_filename):  # Vérifie si le fichier existe déjà
        counter += 1  # Incrémente le compteur
        output_filename = f"{base_name}{counter}.txt"  # Génère un nouveau nom avec le compteur

    return output_filename  # Retourne le nom de fichier disponible

# Fonction pour enregistrer les résultats dans un fichier texte
def save_output_to_file(output_filename, repeats):
    """Enregistre les résultats dans un fichier texte."""
    with open(output_filename, "w") as f:  # Ouvre le fichier en mode écriture
        for repeat in repeats:  # Parcourt chaque répétition maximale trouvée
            word, length, positions, encodings = repeat  # Récupère les informations de la répétition
            # Formatte les positions et les encadrements dans une chaîne de caractères
            enc_pos = ", ".join([f"{pos} ({enc[0]},{enc[1]})" for pos, enc in zip(positions, encodings)])
            # Écrit la répétition, sa longueur, ses positions et ses encadrements dans le fichier
            f.write(f"{word} {length} {enc_pos}\n")

# Fonction principale du programme
def main():
    parser = argparse.ArgumentParser(description="Détection des répétitions maximales dans une séquence ADN.")  # Crée un parser pour gérer les arguments en ligne de commande
    parser.add_argument("fasta_file", help="Fichier FASTA contenant la séquence ADN")  # Argument pour le fichier FASTA
    parser.add_argument("min_len", type=int, nargs='?', default=20, help="Longueur minimale des répétitions (par défaut 20)")  # Argument optionnel pour la longueur minimale des répétitions
    parser.add_argument("max_len", type=int, nargs='?', help="Longueur maximale des répétitions (optionnel)")  # Argument optionnel pour la longueur maximale des répétitions
    parser.add_argument("-c", "--no-overlap", action="store_true", help="Exclure les répétitions chevauchantes")  # Option pour exclure les répétitions chevauchantes

    args = parser.parse_args()  # Analyse les arguments passés en ligne de commande

    # Générer le nom de fichier de sortie
    script_name = os.path.basename(__file__)  # Obtenir le nom du script actuel
    output_filename = generate_output_filename(script_name)  # Génère le nom de fichier de sortie

    # Lecture du fichier FASTA
    seq_id, sequence = parse_fasta(args.fasta_file)  # Lit la séquence ADN du fichier FASTA

    # Recherche des répétitions maximales
    repeats = find_repeats(sequence, args.min_len, args.max_len, args.no_overlap)  # Trouve les répétitions maximales

    # Enregistrer les résultats dans un fichier texte
    save_output_to_file(output_filename, repeats)  # Sauvegarde les répétitions trouvées dans un fichier

    # Indiquer à l'utilisateur où se trouvent les résultats
    print(f"Résultats enregistrés dans {output_filename}")  # Affiche un message avec le nom du fichier de sortie

# Point d'entrée du programme
if __name__ == "__main__":
    main()  # Appelle la fonction principale lorsque le script est exécuté



