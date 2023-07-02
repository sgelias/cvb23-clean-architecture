from pathlib import Path
import click
from decli import cli
from cvb23.core.domain.dtos.blast_config import (
    BlastConfigDTO,
    ReferenceDatabaseDTO,
)

from cvb23.core.domain.dtos.data_type import DataType
from cvb23.core.domain.dtos.query_sequences import QuerySequencesDTO
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
def blast_cmd(
    query: str,
    query_data_type: str,
    database: str,
) -> None:
    try:
        raise NotImplementedError()

        query_path = Path(query)

        query_sequences = QuerySequencesDTO.create(
            path=query_path,
            data_type=DataType(query_data_type),
        )

        reference_database = ReferenceDatabaseDTO.create()

        """ blast_config = BlastConfigDTO.create(

        ) """

        click.echo(query_sequences)

    except Exception as e:
        exc = BlastCmdException(e)
        click.echo(exc, err=True)
        raise exc
