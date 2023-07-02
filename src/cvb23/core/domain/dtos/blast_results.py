from pandera import DataFrameModel, Field


class BlastResultDataTable(DataFrameModel):
    query: str = Field(
        coerce=True,
        nullable=False,
    )

    subject: str = Field(
        coerce=True,
        nullable=False,
    )

    perc_identity: float = Field(
        unique=False,
        coerce=True,
        ge=0.0,
        le=100.0,
        nullable=False,
    )

    align_length: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    mismatches: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    gap_openings: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    q_start: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    q_end: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    s_start: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    s_end: int = Field(
        unique=False,
        coerce=True,
        ge=0,
        nullable=False,
    )

    e_value: float = Field(
        unique=False,
        coerce=True,
        ge=0.0,
        nullable=False,
    )

    bit_score: float = Field(
        unique=False,
        coerce=True,
        ge=0.0,
        nullable=False,
    )
