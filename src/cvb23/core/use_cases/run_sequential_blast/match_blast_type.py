from cvb23.core.domain.dtos.data_type import DataType
from cvb23.core.domain.entities.blast_executor import BlastType
from cvb23.core.domain.exceptions import MatchBlastTypeException


def match_blast_type(
    query_type: DataType,
    subject_type: DataType,
) -> BlastType:
    """Match BLAST type based on query and subject data types.

    Args:
        query_type (DataType): Query data type.
        subject_type (DataType): Subject data type.

    Returns:
        BlastType: BLAST type.

    Raises:
        MatchBlastTypeException: If an error occurs during BLAST type matching.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? 1 - Validate entry parameters
        # ? --------------------------------------------------------------------

        if not isinstance(query_type, DataType):
            raise TypeError("query_type must be of type DataType")

        if not isinstance(subject_type, DataType):
            raise TypeError("subject_type must be of type DataType")

        # ? --------------------------------------------------------------------
        # ? 2 - Execute BLAST
        # ? --------------------------------------------------------------------

        COMBINATIONS = {
            (DataType.NUCLEOTIDE, DataType.NUCLEOTIDE): BlastType.BLASTN,
            (DataType.NUCLEOTIDE, DataType.AMINO_ACID): BlastType.TBLASTX,
            (DataType.AMINO_ACID, DataType.NUCLEOTIDE): BlastType.TBLASTN,
            (DataType.AMINO_ACID, DataType.AMINO_ACID): BlastType.BLASTP,
        }

        if (query_type, subject_type) in COMBINATIONS:
            return COMBINATIONS[(query_type, subject_type)]

        # ? --------------------------------------------------------------------
        # ? 3 - Return a negative response if no match is found
        # ? --------------------------------------------------------------------

        raise MatchBlastTypeException("No Blast type match found")

    except Exception as e:
        raise MatchBlastTypeException(e)
