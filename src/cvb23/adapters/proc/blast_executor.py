from subprocess import check_output

from Bio.SeqRecord import SeqRecord
from pandas import read_csv

from cvb23.core.domain.dtos.blast_config import BlastConfigDTO
from cvb23.core.domain.dtos.blast_results import BlastResultDataTable
from cvb23.core.domain.entities.blast_executor import BlastExecutor, BlastType
from cvb23.core.domain.exceptions import BlastExecutorException


class BlastExecutorProcRepository(BlastExecutor):
    # ? ------------------------------------------------------------------------
    # ? ABSTRACT METHODS IMPLEMENTATIONS
    # ? ------------------------------------------------------------------------

    def run(
        self,
        query_sequences: list[SeqRecord],
        blast_config: BlastConfigDTO,
        blast_type: BlastType,
    ) -> BlastResultDataTable:
        try:
            string_sequences = "\n".join(
                [
                    f">{seq_record.description}\n{seq_record.seq}"
                    for seq_record in query_sequences
                ]
            )

            config = [
                f"{blast_type.value}",
                "-subject",
                f"{blast_config.reference_database.path}",
                "-out",
                f"{blast_config.output_file}",
                "-outfmt",
                f"{blast_config.outfmt}",
            ]

            if blast_config.e_value is not None:
                blast_config.append(f"-evalue {blast_config.e_value}")

            if blast_config.perc_identity is not None:
                blast_config.append(
                    f"-perc_identity {blast_config.perc_identity}"
                )

            if blast_config.max_target_seqs is not None:
                blast_config.append(
                    f"-max_target_seqs {blast_config.max_target_seqs}"
                )

            check_output(
                config,
                input=string_sequences,
                encoding="ascii",
            )

            return BlastResultDataTable.validate(
                read_csv(
                    blast_config.output_file,
                    sep="\t",
                    header=None,
                    names=[
                        BlastResultDataTable.query,
                        BlastResultDataTable.subject,
                        BlastResultDataTable.perc_identity,
                        BlastResultDataTable.align_length,
                        BlastResultDataTable.mismatches,
                        BlastResultDataTable.gap_openings,
                        BlastResultDataTable.q_start,
                        BlastResultDataTable.q_end,
                        BlastResultDataTable.s_start,
                        BlastResultDataTable.s_end,
                        BlastResultDataTable.e_value,
                        BlastResultDataTable.bit_score,
                    ],
                )
            )

        except Exception as e:
            raise BlastExecutorException(e)
