from abc import ABCMeta, abstractmethod
from enum import Enum

from Bio.SeqRecord import SeqRecord

from cvb23.core.domain.dtos.blast_config import BlastConfigDTO
from cvb23.core.domain.dtos.blast_results import BlastResultDataTable


class BlastType(Enum):
    BLASTN = "blastn"
    BLASTP = "blastp"
    BLASTX = "blastx"
    TBLASTN = "tblastn"
    TBLASTX = "tblastx"


class BlastExecutor(metaclass=ABCMeta):
    # ? ------------------------------------------------------------------------
    # ? ABSTRACT METHODS
    # ? ------------------------------------------------------------------------

    @abstractmethod
    def run(
        self,
        query_sequences: list[SeqRecord],
        blast_config: BlastConfigDTO,
        blast_type: BlastType,
    ) -> BlastResultDataTable:
        raise NotImplementedError
