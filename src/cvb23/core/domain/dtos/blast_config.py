from pathlib import Path
from typing import Self

from attrs import define, field

from cvb23.core.domain.dtos.data_type import DataType
from cvb23.core.domain.exceptions import BlastConfigException


@define
class ReferenceDatabaseDTO:
    name: str = field()
    path: Path = field()
    data_type: DataType = field()

    @classmethod
    def create(
        cls,
        path: Path,
        data_type: DataType,
        name: str | None = None,
    ) -> Self:
        try:
            return cls(
                name=path.name if name is None else name,
                path=path,
                data_type=data_type,
            )

        except Exception as e:
            raise BlastConfigException(e)


@define
class BlastConfigDTO:
    # IO Parameters
    reference_database: ReferenceDatabaseDTO = field()
    output_file: Path = field()

    # BLAST Parameters
    outfmt: int = field(default=6)
    e_value: float | None = field(default=None)
    num_threads: int | None = field(default=None)
    word_size: int | None = field(default=None)
    max_target_seqs: int | None = field(default=None)
    perc_identity: float | None = field(default=None)

    @classmethod
    def create(
        cls,
        reference_database: ReferenceDatabaseDTO,
        output_file: Path,
        outfmt: int = 6,
        e_value: float | None = None,
        num_threads: int | None = None,
        word_size: int | None = None,
        max_target_seqs: int | None = None,
        perc_identity: float | None = None,
    ) -> Self:
        try:
            return cls(
                reference_database=reference_database,
                output_file=output_file,
                outfmt=outfmt,
                e_value=e_value,
                num_threads=num_threads,
                word_size=word_size,
                max_target_seqs=max_target_seqs,
                perc_identity=perc_identity,
            )

        except Exception as e:
            raise BlastConfigException(e)
