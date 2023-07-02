from pandera import DataFrameModel, Field


class BlastResultDataTable(DataFrameModel):
    query: str = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        nullable=False,
        strip=True,
    )

    subject: str = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        nullable=False,
        strip=True,
    )

    perc_identity: float = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0.0,
        le=100.0,
        nullable=False,
        strip=True,
    )

    align_length: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    mismatches: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    gap_openings: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    q_start: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    q_end: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    s_start: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    s_end: int = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0,
        nullable=False,
        strip=True,
    )

    e_value: float = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0.0,
        nullable=False,
        strip=True,
    )

    bit_score: float = Field(
        allow_duplicates=False,
        check_dtype=False,
        coerce=True,
        ge=0.0,
        nullable=False,
        strip=True,
    )
