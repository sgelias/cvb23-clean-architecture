from enum import Enum
from os import getenv

from cvb23.core.domain.dtos.data_type import DataType


# ? ----------------------------------------------------------------------------
# ? Get the number of threads from the environment
# ? ----------------------------------------------------------------------------


ENV_NUM_THREADS = getenv("NUM_THREADS", 1)


# ? ----------------------------------------------------------------------------
# ? Configure available databases
# ? ----------------------------------------------------------------------------


class DefaultBlastDatabases(Enum):
    BSUB_GYRB = (
        DataType.AMINO_ACID,
        "src/cvb23/assets/input/ref_dbs/uniprotkb_gene_gyrb_AND_taxonomy_id_492_2023_07_01.faa",
    )
