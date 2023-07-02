from pathlib import Path
from typing import Self

from attrs import define, field

from cvb23.core.domain.dtos.data_type import DataType
from cvb23.core.domain.exceptions import QuerySequenceException


@define
class QuerySequencesDTO:
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
            raise QuerySequenceException(e)
