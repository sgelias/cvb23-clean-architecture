from pathlib import Path

import click
from cvb23.adapters.proc.blast_executor import BlastExecutorProcRepository

from cvb23.core.domain.dtos.blast_config import (
    BlastConfigDTO,
    ReferenceDatabaseDTO,
)
from cvb23.core.domain.dtos.data_type import DataType
from cvb23.core.domain.dtos.query_sequences import QuerySequencesDTO
from cvb23.core.use_cases.run_sequential_blast import run_sequential_blast
from cvb23.ports.cli.exceptions import BlastCmdException
from cvb23.settings import DefaultBlastDatabases


@click.group()
def cvb_cmd() -> None:
    pass


@cvb_cmd.command(
    "blast",
    help="Execute blast through cvb23 wrapper.",
)
@click.option(
    "-q",
    "--query",
    required=True,
    prompt=True,
    type=click.Path(
        resolve_path=False,
        readable=True,
        exists=False,
    ),
    help="The system path of the query file.",
)
@click.option(
    "-qdt",
    "--query-data-type",
    required=True,
    type=click.Choice([i.value for i in DataType]),
    default=DataType.NUCLEOTIDE.value,
    show_default=True,
    help="The query data type.",
)
@click.option(
    "-db",
    "--database",
    required=False,
    type=click.Choice([i.name for i in DefaultBlastDatabases]),
    default=DataType.NUCLEOTIDE.value,
    show_default=True,
    help="Databases to select from.",
)
@click.option(
    "-o",
    "--output-file",
    required=True,
    prompt=True,
    type=click.Path(
        resolve_path=False,
        readable=True,
        exists=False,
    ),
    help="The system path of the output file.",
)
def blast_cmd(
    query: str,
    query_data_type: str,
    database: str,
    output_file: str,
) -> None:
    try:
        # ? --------------------------------------------------------------------
        # ? Build sequences object
        # ? --------------------------------------------------------------------

        query_sequences = QuerySequencesDTO.create(
            path=Path(query),
            data_type=DataType(query_data_type),
        )

        # ? --------------------------------------------------------------------
        # ? Build database object
        # ? --------------------------------------------------------------------

        database_data_type: DataType
        database_path: Path
        (database_data_type, database_path) = next(
            i.value for i in DefaultBlastDatabases if i.name == database
        )

        reference_database = ReferenceDatabaseDTO.create(
            path=database_path,
            data_type=database_data_type,
        )

        # ? --------------------------------------------------------------------
        # ? Build blast configuration object
        # ? --------------------------------------------------------------------

        blast_config = BlastConfigDTO.create(
            reference_database=reference_database,
            output_file=Path(output_file),
        )

        # ? --------------------------------------------------------------------
        # ? Execute blast
        # ? --------------------------------------------------------------------

        blast_response = run_sequential_blast(
            query_sequences=query_sequences,
            blast_config=blast_config,
            blast_executor=BlastExecutorProcRepository(),
        )

        # ? --------------------------------------------------------------------
        # ? Return the user feedback
        # ? --------------------------------------------------------------------

        if blast_response is True:
            click.echo("Blast executed successfully!")
            return

        raise BlastCmdException("Blast executed with errors!")

    except Exception as e:
        exc = BlastCmdException(e)
        click.echo(exc, err=True)
        raise exc
