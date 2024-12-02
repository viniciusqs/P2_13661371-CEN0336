def contar_nos(arvore):
    """
    Função recursiva para contar o número total de nós em uma árvore representada como dicionário de dicionários.

    :param arvore: Dicionário representando a árvore filogenética.
    :return: Número total de nós na árvore.
    """
    total_nos = 1  # Conta o nó atual
    for subarvore in arvore.values():
        total_nos += contar_nos(subarvore)  # Soma os nós das subárvores
    return total_nos


# Testando com a árvore fornecida
arvore = {
    "Filo1": {
        "Classe1": {
            "Ordem1": {
                "Familia1": {},
                "Familia2": {}
            },
            "Ordem2": {
                "Familia3": {
                    "Genero3": {},
                    "Genero4": {}
                }
            }
        },
        "Classe2": {
            "Ordem3": {},
            "Ordem4": {
                "Familia4": {},
                "Familia5": {
                    "Genero1": {},
                    "Genero2": {
                        "Especie1": {},
                        "Especie2": {}
                    }
                }
            }
        }
    }
}

# Chamada da função e impressão do resultado
print("Número total de nós na árvore:", contar_nos(arvore))
