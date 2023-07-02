# CVB23 Blast Wrapper


## The application core resides in src/cvb23/core directory:

```bash
── core
   ├── __init__.py
   ├── domain
   │   ├── __init__.py
   │   ├── dtos
   │   │   ├── __init__.py
   │   │   ├── blast_config.py
   │   │   └── blast_results.py
   │   ├── entities
   │   │   ├── __init__.py
   │   │   └── run_solo_blast.py
   │   └── exceptions
   │       └── __init__.py
   └── use_cases
       ├── __init__.py
       ├── build_consensus_from_tabular_results
       │   └── __init__.py
       └── run_sequential_blast
           └── __init__.py
```
