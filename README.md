# CVB23 Blast Wrapper

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [CVB23 Blast Wrapper](#cvb23-blast-wrapper)
  - [Estrutura do projeto](#estrutura-do-projeto)
  - [Utilização da ferramenta](#utilização-da-ferramenta)
    - [Instalação](#instalação)
    - [Execução](#execução)

<!-- /code_chunk_output -->

## Estrutura do projeto

Esse é o projeto base utilizado no CVB23. O projeto é dividido em três partes
principais:

1. Core (src/cvb23/core):
   * Aqui estão as regras de negócio do projeto. Elas incluem os principais
      casos de uso, os objetos de trânsito de informação (DTOs) e as entidades
      do projeto.

```bash
src/cvb23/core
├── domain
│   ├── dtos
│   │   ├── blast_config.py
│   │   ├── blast_results.py
│   │   ├── data_type.py
│   │   ├── __init__.py
│   │   └── query_sequences.py
│   ├── entities
│   │   ├── blast_executor.py
│   │   └── __init__.py
│   ├── exceptions
│   │   └── __init__.py
│   └── __init__.py
├── __init__.py
└── use_cases
    ├── __init__.py
    └── run_sequential_blast
        ├── __init__.py
        └── match_blast_type.py
```

1. Adapters (src/cvb23/adapters/proc):
   * Aqui está o adaptador de processo utilizado na execução do blast no
      sistema.

```bash
src/cvb23/adapters
├── __init__.py
└── proc
    ├── blast_executor.py
    └── __init__.py
```

3. Ports (src/cvb23/ports/cli):
   * Aqui está a porta de CLI utilizada para executar o programa.

```bash
src/cvb23/ports
├── cli
│   ├── exceptions.py
│   ├── __init__.py
│   └── main.py
└── __init__.py
```

## Utilização da ferramenta

### Instalação

Para executar a ferramenta é necessário realizar sua instalação prévia
utilizando a ferramenta pipenv (para testes) ou diretamente seu ambiente nativo
com o interpretador com suporte mínimo de python3.11.

Para instalar utilizando Pipenv, execute o comando abaixo no terminal Linux,
estando dentro da pasta raiz do projeto:

```bash
pipenv install --dev -e '.'
```

Para instalar diretamente no seu ambiente utilize o comando abaixo:

```bash
python3.11 -m pip install -e .
```


### Execução


```bash
cvblast-cmd blast \
   -q src/cvb23/assets/input/query/query.faa \
   -qdt amino_acid \
   -db BSUB_GYRB \
   -o src/cvb23/assets/output/blast.out
```
