import sys
import os

def find_longest_orf(seq, frame):
    """
    Encontra o ORF mais longo em uma sequência de DNA para o quadro especificado.

    Args:
        seq (str): Sequência de DNA.
        frame (int): Quadro de leitura inicial (0, 1 ou 2).

    Returns:
        str: O ORF mais longo encontrado.
    """
    # Códons que iniciam e terminam uma fase de leitura
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    longest_orf = ""  # Armazena o maior ORF encontrado
    start = -1  # Posição inicial do ORF; -1 indica que ainda não iniciou

    # Percorre a sequência em intervalos de 3 (códons), começando no quadro especificado
    for i in range(frame, len(seq), 3):
        codon = seq[i:i+3]  # Extrai o códon atual (3 bases)
        if len(codon) < 3:  # Interrompe se o códon estiver incompleto no final
            break
        if codon == start_codon and start == -1:
            start = i  # Marca o início do ORF
        elif codon in stop_codons and start != -1:
            # ORF encontrado entre o códon de início e o de parada
            orf = seq[start:i+3]
            if len(orf) > len(longest_orf):  # Atualiza se for o maior encontrado
                longest_orf = orf
            start = -1  # Reinicia a busca por um novo ORF
    return longest_orf

def translate_orf(orf):
    """
    Traduza uma sequência de DNA para um peptídeo usando uma tabela de códons.

    Args:
        orf (str): Sequência de DNA representando um ORF.

    Returns:
        str: Sequência de aminoácidos (peptídeo).
    """
    # Tabela de códons: associa tríades de nucleotídeos a aminoácidos
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W'
    }
    peptide = ""  # Armazena a sequência de aminoácidos gerada
    # Percorre o ORF em grupos de 3 bases (códons)
    for i in range(0, len(orf), 3):
        codon = orf[i:i+3]  # Extrai o códon atual
        if len(codon) == 3:  # Tradução só ocorre para códons completos
            peptide += codon_table.get(codon, 'X')  # 'X' para códons desconhecidos
    return peptide

def parse_fasta(input_file):
    """
    Lê um arquivo FASTA e separa os registros.

    Args:
        input_file (str): Caminho para o arquivo FASTA.

    Returns:
        list: Lista de tuplas (header, sequence).
    """
    fasta_records = []  # Lista para armazenar os registros (header, sequence)
    try:
        # Lê o arquivo e separa os registros pelo caractere '>'
        with open(input_file, "r") as file:
            data = file.read().strip().split(">")[1:]  # Ignora tudo antes do primeiro '>'
            for record in data:
                lines = record.splitlines()
                header = lines[0]  # A primeira linha é o identificador da sequência
                sequence = ''.join(lines[1:])  # Junta todas as linhas da sequência
                fasta_records.append((header, sequence))
    except Exception as e:
        raise ValueError(f"Erro ao processar o arquivo FASTA: {e}")
    return fasta_records

def process_records(records, output_orf, output_peptide):
    """
    Processa os registros FASTA, encontra ORFs e escreve os resultados.

    Args:
        records (list): Lista de tuplas (header, sequence).
        output_orf (str): Caminho para o arquivo de saída dos ORFs.
        output_peptide (str): Caminho para o arquivo de saída dos peptídeos.
    """
    # Abre os arquivos de saída para escrita
    with open(output_orf, "w") as orf_out, open(output_peptide, "w") as pep_out:
        for header, sequence in records:
            longest_orf = ""  # Armazena o maior ORF encontrado
            longest_frame = 0  # Quadro de leitura correspondente ao maior ORF
            start_pos = 0  # Posição inicial do maior ORF
            end_pos = 0  # Posição final do maior ORF

            # Busca o maior ORF em cada um dos três quadros de leitura
            for frame in range(3):
                orf = find_longest_orf(sequence, frame)
                if len(orf) > len(longest_orf):
                    longest_orf = orf
                    longest_frame = frame + 1
                    start_pos = sequence.find(orf) + 1
                    end_pos = start_pos + len(orf) - 1

            if longest_orf:
                # Cria o identificador para o maior ORF
                output_header = f"{header}_frame{longest_frame}_{start_pos}_{end_pos}"
                # Escreve o ORF no arquivo de saída
                orf_out.write(f">{output_header}\n{longest_orf}\n")
                # Tradução do ORF para peptídeo
                peptide = translate_orf(longest_orf)
                # Escreve o peptídeo no arquivo de saída
                pep_out.write(f">{output_header}\n{peptide}\n")

if __name__ == "__main__":
    try:
        # Define o arquivo de entrada e os arquivos de saída
        input_file = sys.argv[1] if len(sys.argv) > 1 else "input_example.fasta"
        output_orf_file = "ORF.fna"
        output_peptide_file = "ORF.faa"

        # Verifica se o arquivo de entrada existe
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Arquivo {input_file} não encontrado.")

        # Processa o arquivo FASTA
        fasta_records = parse_fasta(input_file)

        # Processa os registros e escreve os resultados
        process_records(fasta_records, output_orf_file, output_peptide_file)

        print(f"Processamento concluído. Arquivos gerados: {output_orf_file}, {output_peptide_file}")

    except Exception as e:
        # Captura qualquer erro inesperado durante a execução
        print(f"Erro durante a execução: {e}")
