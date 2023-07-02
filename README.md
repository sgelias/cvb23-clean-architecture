# CVB23 Blast Wrapper

## Descrição

Esse é o projeto base utilizado no CVB23. O projeto é dividido em três partes
principais:

1. Core (src/cvb23/core):
   * Aqui estão as regras de negócio do projeto. Elas incluem os principais
      casos de uso, os objetos de trânsito de informação (DTOs) e as entidades
      do projeto.

## O core da aplicação src/cvb23/core directory:

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

2. Adapters (src/cvb23/adapters/proc):
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
