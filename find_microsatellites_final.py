# microsatellite_info_acgt_updated_graph.py
import re
import matplotlib.pyplot as plt

def find_acgt_microsatellites(sequence):
    acgt_microsatellites = re.finditer(r'(ACGT){3,}', sequence)
    positions = [(match.start(), len(match.group())) for match in acgt_microsatellites]
    return positions

def count_acgt_repeats(sequence):
    return len(re.findall(r'(ACGT){3,}', sequence))

def count_acgt_motifs(acgt_microsatellites):
    motif_counts = {}
    for position, length in acgt_microsatellites:
        motif_counts[(position, length)] = motif_counts.get((position, length), 0) + 1
    return motif_counts

def plot_acgt_microsatellite_info(acgt_microsatellites, acgt_repeat_count):
    positions = [position for position, length in acgt_microsatellites]
    counts = [motif_count for _, motif_count in acgt_microsatellites]

    y_labels = [f'Position: {position}\nRepeats: {count}' for position, count in zip(positions, counts)]

    plt.barh(y_labels, counts)
    plt.xlabel('Count of "ACGT" Motif Repeats')
    plt.ylabel('Position of "ACGT" Microsatellite')
    plt.title('"ACGT" Microsatellite Position and Repeats Distribution')

    plt.show()

def read_sequence_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            sequence = file.read().replace('\n', '')
            if not sequence:
                print("Empty genetic sequence")
            return sequence
    except FileNotFoundError:
        print("File not found.")
        return None

def main():
    file_path = input("Enter the path to the file containing the genetic sequence: ")
    sequence = read_sequence_from_file(file_path)

    if sequence is not None:
        acgt_microsatellites = find_acgt_microsatellites(sequence)
        acgt_repeat_count = count_acgt_repeats(sequence)

        if acgt_microsatellites:
            print('\n"ACGT" Microsatellite SSR Motifs:')
            for position, length in acgt_microsatellites:
                print(f'Position: {position}, Length: {length}')

            motif_counts = count_acgt_motifs(acgt_microsatellites)
            print('\nCount of "ACGT" Motifs for Each Microsatellite:')
            for (position, length), count in motif_counts.items():
                print(f'Position: {position}, Length: {length}, Count: {count}')

            print(f'\nNumber of "ACGT" Microsatellite Repeats: {acgt_repeat_count}')
            plot_acgt_microsatellite_info(list(motif_counts.keys()), acgt_repeat_count)
        else:
            print("Sequence with no (ACGT)n microsatellites")
    else:
        print("File with no genetic sequence")

if __name__ == "__main__":
    main()
