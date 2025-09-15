import logging
from collections import defaultdict
import json
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import TableConfig


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module parses output from the [Space Ranger pipeline]("https://www.10xgenomics.com/support/software/space-ranger/latest").
    Space Ranger is a set of pipelines used to process 10X Visium data. This module parses both the metrics_summary.csv and web_summary.html
    files output from the Space Ranger pipeline in order to recover QC metrics about a given analysis. This module supports output from Space Ranger
    1.0.0 through 4.0.1.

    Output from Space Ranger typically looks like this:
    `
    Outputs:
    - Run summary HTML:                 /home/jdoe/runs/sample345/outs/web_summary.html
    - Outputs of spatial pipeline:
                aligned_fiducials:              /home/jdoe/runs/sample345/outs/spatial/aligned_fiducials.jpg
                detected_tissue_image:          /home/jdoe/runs/sample345/outs/spatial/detected_tissue_image.jpg
                scalefactors_json:              /home/jdoe/runs/sample345/outs/spatial/scalefactors_json.json
                tissue_hires_image:             /home/jdoe/runs/sample345/outs/spatial/tissue_hires_image.png
                tissue_lowres_image:            /home/jdoe/runs/sample345/outs/spatial/tissue_lowres_image.png
                cytassist_image:                null
                aligned_tissue_image:           null
                tissue_positions:               /home/jdoe/runs/sample345/outs/spatial/tissue_positions.csv
                spatial_enrichment:             /home/jdoe/runs/sample345/outs/spatial/spatial_enrichment.csv
                barcode_fluorescence_intensity: null
        - Run summary CSV:                                           /home/jdoe/runs/sample345/outs/metrics_summary.csv
        - Correlation values between isotypes and Antibody features: null
        - BAM:                                                       /home/jdoe/runs/sample345/outs/possorted_genome_bam.bam
        - BAM BAI index:                                             /home/jdoe/runs/sample345/outs/possorted_genome_bam.bam.bai
        - BAM CSI index:                                             null
        - Filtered feature-barcode matrices MEX:                     /home/jdoe/runs/sample345/outs/filtered_feature_bc_matrix
        - Filtered feature-barcode matrices HDF5:                    /home/jdoe/runs/sample345/outs/filtered_feature_bc_matrix.h5
        - Unfiltered feature-barcode matrices MEX:                   /home/jdoe/runs/sample345/outs/raw_feature_bc_matrix
        - Unfiltered feature-barcode matrices HDF5:                  /home/jdoe/runs/sample345/outs/raw_feature_bc_matrix.h5
        - Secondary analysis output CSV:                             /home/jdoe/runs/sample345/outs/analysis
        - Per-molecule read information:                             /home/jdoe/runs/sample345/outs/molecule_info.h5
        - Loupe Browser file:                                        /home/jdoe/runs/sample345/outs/cloupe.cloupe
        - Feature Reference:                                         null
        - Target Panel file:                                         null
        - Probe Set file:                                            null
        Pipestance completed successfully!
    `
    Only a metrics summary file OR a web summary file is required. Both can be provided to recover more information about an analysis,
    as these files potentially report different metrics (depending on the version of Space Ranger/assay used). For example, in a
    Visium HD analysis metrics about 2µm bins are only provided in the metrics_summary.csv. This information is not available in
    the web summary.

    h/t: https://github.com/MultiQC/MultiQC/blob/main/multiqc/modules/xenium/xenium.py
    *Todo? Todone?
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Space Ranger",
            anchor="spaceranger",
            href=["https://www.10xgenomics.com/support/software/space-ranger/latest"],
            info="Space Ranger is a set of analysis pipelines that process 10x Genomics Visium data with brightfield or fluorescence microscope images, allowing users to map the whole transcriptome in a variety of tissues",
        )

        data_by_sample = {}
        for f in self.find_log_files("spaceranger/metrics", filehandles=True):
            parsed_data = self.parse_spaceranger_metrics(f)
            if parsed_data:
                sample_name = parsed_data["Sample ID"]
                if sample_name in data_by_sample.keys():
                    log.critical(
                        f"Sample name '{sample_name}' is shared across multiple metrics summary files in this run. Overwriting {sample_name}."
                    )
                data_by_sample[sample_name] = parsed_data
                self.add_data_source(f, sample_name)

        data_by_sample_html = {}
        for f2 in self.find_log_files("spaceranger/count_html", filehandles=True):
            parsed_data_count = self.parse_count_html(f2)
            if parsed_data_count:
                sample_name = parsed_data_count["Sample ID"]
                if sample_name in data_by_sample_html.keys():
                    log.critical(
                        f"Sample name '{sample_name}' is shared across multiple count html files in this run. Overwriting {sample_name}."
                    )
                data_by_sample_html[sample_name] = parsed_data_count
                try:
                    # This should have the CSV values overwrite the ones from the HTML file
                    data_by_sample[sample_name] = {**data_by_sample_html[sample_name], **data_by_sample[sample_name]}
                except KeyError:
                    log.warning(
                        f"Sample {sample_name} does not seem to have a metrics summary file, pulling data from the HTML file instead."
                    )
                    data_by_sample[sample_name] = data_by_sample_html[sample_name]
                self.add_data_source(f2, sample_name)

        data_by_sample = self.ignore_samples(data_by_sample)
        log.info(f"Found {len(data_by_sample)} Space Ranger reports")
        self.write_data_file(data_by_sample, "multiqc_spaceranger")
        self.spaceranger_general_stats_table(data_by_sample)

        self.add_section(
            name="Analysis results",
            anchor="Table",
            description="Summary Stats from Space Ranger run",
            helptext="""
                    Summary Stats from Space Ranger run
                    """,
            plot=self.spaceranger_transcript_table(data_by_sample),
        )

        self.add_section(
            name="Sequencing Saturation",
            anchor="Sequencing Saturation",
            description="Plot of sequencing saturation",
            helptext="""
                    Plot of sequencing saturation
                    """,
            plot=self.add_seq_sat_plot(data_by_sample),
        )

        self.add_section(
            name="Genomic UMIs",
            anchor="Genomic UMIs",
            description="Plot of Genomic UMIs",
            helptext="""
                    Plot of Genomic UMIs
                    """,
            plot=self.add_gen_umi_plot(data_by_sample),
        )

        self.add_section(
            name="Genes detected",
            anchor="Plot",
            description="Plot of Genes detected in differing sized bins",
            helptext="""
                    Plot of Genes detected in differing sized bins
                    """,
            plot=self.add_gene_number_plot(data_by_sample),
        )

    def add_seq_sat_plot(self, data_by_sample):
        config = {"ylab": "Sequencing Saturation (%)", "cpswitch": False}
        config["id"] = "Space Ranger Sequencing Saturation plot"
        config["title"] = "Space Ranger: Saturation plot"

        genes_detected = {"Sequencing Saturation": {"color": "#f7a35c", "name": "Sequencing Saturation"}}
        return bargraph.plot(data_by_sample, genes_detected, config)

    def add_gen_umi_plot(self, data_by_sample):
        config = {"ylab": "Fraction Genomic UMIs (%)", "cpswitch": False}

        config["id"] = "Space Ranger Genomic UMIs plot"
        config["title"] = "Space Ranger: Genomic UMIs"

        genes_detected = {
            "Estimated UMIs from Genomic DNA": {"color": "#320faf", "name": "Estimated UMIs from Genomic DNA"}
        }
        return bargraph.plot(data_by_sample, genes_detected, config)

    def add_gene_number_plot(self, data_by_sample):
        config = {"ylab": "Genes detected by bin size", "stacking": "group", "cpswitch": False}
        config["id"] = "Space Ranger: Genes detected by bin size"
        config["title"] = "Space Ranger: Genes detected by bin size"
        genes_detected = {
            "Mean Genes Under Tissue per Square 2 µm": {
                "color": "#20568f",
                "name": "Mean Genes Under Tissue per Square 2 µm",
            },
            "Mean Genes Under Tissue per Bin 8 µm": {
                "color": "#f7a35c",
                "name": "Mean Genes Under Tissue per Bin 8 µm",
            },
            "Mean Genes Under Tissue per Bin 16 µm": {
                "color": "#981919",
                "name": "Mean Genes Under Tissue per Bin 16 µm",
            },
        }

        return bargraph.plot(data_by_sample, genes_detected, config)

    def parse_spaceranger_metrics(self, f):
        """
        Parse Space Ranger metrics_summary.csv file

        Parses sample metrics from metrics_summary.csv file. Notably, for Visium HD this file contains 2µm metrics that are
        not present in the count HTML file. Also much easier to parse and probably won't change if 10X decide to restructure
        their web reports. Probably.

        """
        lines = f["f"].read().splitlines()
        if len(lines) < 2:
            return {}

        # Get header and data row
        header = lines[0].split(",")
        data_row = lines[1].split(",")

        if len(header) != len(data_row):
            log.warning(f"Header and data row lengths don't match in {f['fn']}")
            return {}

        # Create mapping of headers to values
        metrics = dict(zip(header, data_row))

        # Convert numeric fields - thoughts on moving this out of a function? it's copied verbatim in the count html parser.
        numeric_fields = [
            # visium hd
            "Number of Reads",
            "Valid Barcodes",
            "Valid UMI Sequences",
            "Sequencing Saturation",
            "Q30 Bases in Barcode",
            "Q30 Bases in Probe Read",
            "Q30 Bases in UMI",
            "Reads Mapped to Probe Set",
            "Reads Mapped Confidently to Probe Set",
            "Fraction Reads in Squares Under Tissue",
            "Genes Detected",
            "Reads Mapped Confidently to the Filtered Probe Set",
            "Number of Genes",
            "Reads Half-Mapped to Probe Set",
            "Reads Split-Mapped to Probe Set",
            "Estimated UMIs from Genomic DNA",
            "Estimated UMIs from Genomic DNA per Unspliced Probe",
            "Number of Squares Under Tissue 2 µm",
            "Mean Reads Under Tissue per Square 2 µm",
            "Fraction of Squares Under Tissue 2 µm",
            "Mean Genes Under Tissue per Square 2 µm",
            "Mean UMIs Under Tissue per Square 2 µm",
            "Total Genes Detected Under Tissue 2 µm",
            "Number of Bins Under Tissue 8 µm",
            "Mean Reads Under Tissue per Bin 8 µm",
            "Fraction of Bins Under Tissue 8 µm",
            "Mean Genes Under Tissue per Bin 8 µm",
            "Mean UMIs Under Tissue per Bin 8 µm",
            "Total Genes Detected Under Tissue 8 µm",
            "UMIs per sq mm of Tissue",
            "Reads per sq mm of Tissue",
            "Number of Bins Under Tissue 16 µm",
            "Mean Reads Under Tissue per Bin 16 µm",
            "Fraction of Bins Under Tissue 16 µm",
            "Mean Genes Under Tissue per Bin 16 µm",
            "Mean UMIs Under Tissue per Bin 16 µm",
            "Total Genes Detected Under Tissue 16 µm",
            "Number of Cells",
            "Reads in Cells",
            "UMIs in Cells",
            "Mean Reads per Cell",
            "Median Genes per Cell",
            "Median UMIs per Cell",
            "Median Cell Area (μm²)Median Nucleus Area (μm²)",
            "Maximum Nucleus Diameter (pixels)",
            # visium SD
            "Number of Spots Under Tissue",
            "Mean Reads per Spot",
            "Mean Reads Under Tissue per Spot",
            "Fraction of Spots Under Tissue",
            "Median Genes per Spot",
            "Median UMI Counts per Spot",
            "Valid Barcodes",
            "Valid UMIs",
            "Sequencing Saturation",
            "Q30 Bases in Barcode",
            "Q30 Bases in RNA Read",
            "Q30 Bases in UMI",
            "Reads Mapped to Genome",
            "Reads Mapped Confidently to Genome",
            "Reads Mapped Confidently to Intergenic Regions",
            "Reads Mapped Confidently to Intronic Regions",
            "Reads Mapped Confidently to Exonic Regions",
            "Reads Mapped Confidently to Transcriptome",
            "Reads Mapped Antisense to Gene",
            "Fraction Reads in Spots Under Tissue",
            "Total Genes Detected",
        ]

        parsed_metrics = {}
        for field in numeric_fields:
            if field in metrics and metrics[field] != "":
                try:
                    parsed_metrics[field] = float(metrics[field])
                except ValueError:
                    log.warning(f"Could not convert {field}='{metrics[field]}' to float")

        # Keep string fields
        string_fields = ["Sample ID"]
        for field in string_fields:
            if field in metrics:
                parsed_metrics[field] = metrics[field]

        # covert some field to percent
        covert_to_percent = ["Sequencing Saturation", "Estimated UMIs from Genomic DNA"]

        if field in parsed_metrics and parsed_metrics[field] != "":
            for field in covert_to_percent:
                if field in parsed_metrics:
                    parsed_metrics[field] = parsed_metrics[field] * 100

        return parsed_metrics

    def parse_count_html(self, f):
        """
        Space Ranger count report parser

        Used to parse additional metrics from Space Ranger count HTML reports.
        The count reports contain some information that cannot be recovered from the metrics summary CSV, such as transcriptome version and
        the pipeline version. If the CSV is absent, all metrics will be pulled from the HTML file instead. Likely to break if 10X decide to
        redesign their web reports. Currently tested with Space Ranger 1.0.0 to 4.0.1.

        This function is adapted from the spaceranger module present in the multiqc 1.31 release.
        """

        warnings_data_by_sample: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict(lambda: defaultdict())

        warnings_headers: Dict = dict()
        summary = None

        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                summary = json.loads(line)
                if "summary" in summary.keys():
                    summary = summary["summary"]
                break

        if summary is None:
            logging.error(f"Couldn't find JSON summary data in HTML report, skipping: {f['fn']}")

        sample_name = self.clean_s_name(summary["sample"]["id"], f)

        # List of data collated from different tables in cellranger reports.
        # This is a list of Tuples (metric name, value)

        # For spacreanger 1 & 2:
        if "summary_tab" in summary.keys():
            software = next(
                iter(x[1] for x in summary["summary_tab"]["pipeline_info_table"]["rows"] if x[0] == "Pipeline Version")
            )

            data_rows = (
                summary["summary_tab"]["pipeline_info_table"]["rows"]
                + summary["summary_tab"]["cells"]["table"]["rows"]
                + summary["summary_tab"]["sequencing"]["table"]["rows"]
                + summary["summary_tab"]["mapping"]["table"]["rows"]
            )

        # For spaceranger 3 & 4:
        if "tabs" in summary.keys():
            software = next(
                iter(
                    x[1]
                    for x in summary["tabs"]["tab_data"][0]["run_summary"]["card"]["inner"]["rows"]
                    if x[0] == "Pipeline Version"
                )
            )

            data_rows = (
                summary["tabs"]["tab_data"][0]["run_summary"]["card"]["inner"]["rows"]
                + summary["tabs"]["tab_data"][3]["segmentation_metrics_table"]["card"]["inner"]["rows"]
                + summary["tabs"]["tab_data"][2]["bin_metrics"]["card"]["inner"]["rows"]
                + summary["tabs"]["tab_data"][0]["top"]["left"]["sequencing_metrics"]["inner"]["rows"]
                + summary["tabs"]["tab_data"][0]["top"]["left"]["mapping_metrics"]["inner"]["rows"]
            )

        # Add software version to report.
        software_name, software_version = software.split("-")
        self.add_software_version(version=software_version, sample=sample_name, software_name=software_name)

        # Extract warnings if any
        alarms_list = summary["alarms"].get("alarms", [])
        for alarm in alarms_list:
            # "Intron mode used" alarm added in Space Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            warnings_data_by_sample[sample_name][alarm["id"]] = "FAIL"
            warnings_headers[alarm["id"]] = {
                "title": alarm["id"].replace("_", " ").title(),
                "description": alarm["title"],
                "bgcols": {"FAIL": "#f7dddc"},
            }

        if len(warnings_data_by_sample) > 0:
            self.add_section(
                name="Count - Warnings",
                anchor="spaceranger-count-warnings-section",
                description="Warnings encountered during the analysis",
                plot=table.plot(
                    warnings_data_by_sample,
                    warnings_headers,
                    {
                        "namespace": "Space Ranger Count",
                        "id": "spaceranger-count-warnings",
                        "title": "Space Ranger: Count: Warnings",
                    },
                ),
            )

        # Convert list of tuples to dict to match the csv metrics
        html_dict = dict(zip([item[0] for item in data_rows], [item[1] for item in data_rows]))

        numeric_fields = [
            # visium hd
            "Number of Reads",
            "Valid Barcodes",
            "Valid UMI Sequences",
            "Sequencing Saturation",
            "Q30 Bases in Barcode",
            "Q30 Bases in Probe Read",
            "Q30 Bases in UMI",
            "Reads Mapped to Probe Set",
            "Reads Mapped Confidently to Probe Set",
            "Fraction Reads in Squares Under Tissue",
            "Genes Detected",
            "Reads Mapped Confidently to the Filtered Probe Set",
            "Number of Genes",
            "Reads Half-Mapped to Probe Set",
            "Reads Split-Mapped to Probe Set",
            "Estimated UMIs from Genomic DNA",
            "Estimated UMIs from Genomic DNA per Unspliced Probe",
            "Number of Squares Under Tissue 2 µm",
            "Mean Reads Under Tissue per Square 2 µm",
            "Fraction of Squares Under Tissue 2 µm",
            "Mean Genes Under Tissue per Square 2 µm",
            "Mean UMIs Under Tissue per Square 2 µm",
            "Total Genes Detected Under Tissue 2 µm",
            "Number of Bins Under Tissue 8 µm",
            "Mean Reads Under Tissue per Bin 8 µm",
            "Fraction of Bins Under Tissue 8 µm",
            "Mean Genes Under Tissue per Bin 8 µm",
            "Mean UMIs Under Tissue per Bin 8 µm",
            "Total Genes Detected Under Tissue 8 µm",
            "UMIs per sq mm of Tissue",
            "Reads per sq mm of Tissue",
            "Number of Bins Under Tissue 16 µm",
            "Mean Reads Under Tissue per Bin 16 µm",
            "Fraction of Bins Under Tissue 16 µm",
            "Mean Genes Under Tissue per Bin 16 µm",
            "Mean UMIs Under Tissue per Bin 16 µm",
            "Total Genes Detected Under Tissue 16 µm",
            "Number of Cells",
            "Reads in Cells",
            "UMIs in Cells",
            "Mean Reads per Cell",
            "Median Genes per Cell",
            "Median UMIs per Cell",
            "Median Cell Area (μm²)Median Nucleus Area (μm²)",
            "Maximum Nucleus Diameter (pixels)",
            # visium SD
            "Number of Spots Under Tissue",
            "Number of Reads",
            "Mean Reads per Spot",
            "Mean Reads Under Tissue per Spot",
            "Fraction of Spots Under Tissue",
            "Median Genes per Spot",
            "Median UMI Counts per Spot",
            "Valid Barcodes",
            "Valid UMIs",
            "Sequencing Saturation",
            "Q30 Bases in Barcode",
            "Q30 Bases in RNA Read",
            "Q30 Bases in UMI",
            "Reads Mapped to Genome",
            "Reads Mapped Confidently to Genome",
            "Reads Mapped Confidently to Intergenic Regions",
            "Reads Mapped Confidently to Intronic Regions",
            "Reads Mapped Confidently to Exonic Regions",
            "Reads Mapped Confidently to Transcriptome",
            "Reads Mapped Antisense to Gene",
            "Fraction Reads in Spots Under Tissue",
            "Total Genes Detected",
        ]

        # Keep string fields
        string_fields = ["Sample ID", "Transcriptome", "Probe Set Name", "Slide Serial Number"]

        parsed_metrics = {}
        for field in numeric_fields:
            # Check for str prior to float conversion, stops excess error messages in log
            if field in html_dict.keys() and isinstance(html_dict[field], str):
                try:
                    # Remove commas and % from numbers in the HTML report
                    parsed_metrics[field] = html_dict[field].replace(",", "").replace("%", "")
                    parsed_metrics[field] = float(parsed_metrics[field])
                except ValueError:
                    log.warning(f"Could not convert {field}='{html_dict[field]}' to float")

                # Correct format of HTML metrics to work with plots
                if (
                    field.startswith("Reads Mapped")
                    or field.startswith("Valid ")
                    or field.startswith("Q30 ")
                    or field.startswith("Fraction Reads")
                ):
                    parsed_metrics[field] = parsed_metrics[field] / 100

        for field in string_fields:
            if field in html_dict.keys():
                parsed_metrics[field] = html_dict[field]

        return parsed_metrics

    def spaceranger_general_stats_table(self, data_by_sample):
        """
        Add key Visium metrics to the general statistics table
        """
        headers: Dict[str, Dict[str, Any]] = {
            "Number of Reads": {
                "title": "Number of Reads",
                "description": "Total number of reads sequenced",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "Chemistry": {
                "title": "Chemistry",
                "description": "Chemistry",
            },
            "Slide Serial Number": {
                "title": "Slide Serial Number",
                "description": "Slide Serial Number",
            },
            "Transcriptome": {
                "title": "Transcriptome",
                "description": "Transcriptome",
            },
            "Probe Set Name": {
                "title": "Probe Set Name",
                "description": "Probe Set Name",
            },
            "Reads Mapped to Probe Set": {
                "title": "Reads Mapped to Probe Set",
                "description": "Reads Mapped to Probe Set",
                "modify": lambda x: x * 100.0,
                "scale": "Purples",
                "format": "{:.2f}%",
            },
            "Reads Mapped to Genome": {
                "title": "Reads Mapped to Genome",
                "description": "Reads Mapped to Genome",
                "modify": lambda x: x * 100.0,
                "scale": "Purples",
                "format": "{:.2f}%",
            },
            "Number of Spots Under Tissue": {
                "title": "Number of Spots Under Tissue",
                "description": "Number of Spots Under Tissue",
                "scale": "Greens",
                "format": "{:,.2f}",
            },
            "Median Genes per Spot": {
                "title": "Median Genes per Spot",
                "description": "Median Genes per Spot",
                "scale": "Blues",
                "format": "{:,.2f}",
            },
            "Reads in Squares Under Tissue": {
                "title": "Reads in Squares Under Tissue",
                "description": "Reads in Squares Under Tissue",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "Mean Genes Under Tissue per Square 2 µm": {
                "title": "Mean Genes Under Tissue per Square 2 µm",
                "description": "Mean Genes Under Tissue per Square 2 µm",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Mean Genes Under Tissue per Bin 8 µm": {
                "title": "Mean Genes Under Tissue per Bin 8 µm",
                "description": "Mean Genes Under Tissue per Bin 8 µm",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "Mean Genes Under Tissue per Bin 16 µm": {
                "title": "Mean Genes Under Tissue per Bin 16 µm",
                "description": "Mean Genes Under Tissue per Bin 16 µm",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Median Genes per Cell": {
                "title": "Median Genes per Cell",
                "description": "Median Genes per Cell",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
        }

        self.general_stats_addcols(data_by_sample, headers)

    def spaceranger_transcript_table(self, data_by_sample):
        headers: Dict[str, Dict[str, Any]] = {
            "Number of Reads": {
                "title": "Number of Reads",
                "description": "Total number of reads sequenced",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Sequencing Saturation": {
                "title": "Sequencing Saturation",
                "description": "Sequencing Saturation",
                "scale": "",
                "format": "{:.2f}%",
            },
            "Valid Barcodes": {
                "title": "Valid Barcodes",
                "description": "Valid Barcodes",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },
            "Valid UMIs": {
                "title": "Valid UMIs",
                "description": "Valid UMIs",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },
            "Reads Mapped to Probe Set": {
                "title": "Reads Mapped to Probe Set",
                "description": "Reads Mapped to Probe Set",
                "modify": lambda x: x * 100.0,
                "scale": "",
                "format": "{:.2f}%",
            },
            "Reads Mapped Confidently to Probe Set": {
                "title": "Reads Mapped Confidently to Probe Set",
                "description": "Reads Mapped Confidently to Probe Set",
                "modify": lambda x: x * 100.0,
                "scale": "",
                "format": "{:.2f}%",
            },
            "Reads Mapped Confidently to the Filtered Probe Set": {
                "title": "Reads Mapped Confidently to the Filtered Probe Set",
                "description": "Reads Mapped Confidently to the Filtered Probe Set",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },
            "Reads Mapped to Genome": {
                "title": "Reads Mapped to Genome",
                "description": "Reads Mapped to Genome",
                "modify": lambda x: x * 100.0,
                "scale": "",
                "format": "{:.2f}%",
            },
            "Reads Mapped Confidently to Genome": {
                "title": "Reads Mapped Confidently to Genome",
                "description": "Reads Mapped Confidently to Genome",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },
            "Reads Mapped Confidently to Exonic Regions": {
                "title": "Reads Mapped Confidently to Exonic Regions",
                "description": "Reads Mapped Confidently to Exonic Regions",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },
            "Total Genes Detected": {
                "title": "Genes Detected (Genome)",
                "description": "Genes Detected (Genome)",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Genes Detected": {
                "title": "Genes Detected (Probe Set)",
                "description": "Genes Detected (Probe Set)",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Median Genes per Spot": {
                "title": "Median Genes per Spot",
                "description": "Median Genes per Spot",
                "scale": "",
                "format": "{:,.2f}",
            },
            "Q30 Bases in Barcode": {
                "title": "Q30 Bases in Barcode",
                "description": "Q30 Bases in Barcode",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:,.2f}%",
            },
            "Q30 Bases in UMI": {
                "title": "Q30 Bases in UMI",
                "description": "Q30 Bases in UMI",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:,.2f}%",
            },
            "Q30 Bases in RNA Read": {
                "title": "Q30 Bases in RNA Read",
                "description": "Q30 Bases in RNA Read",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:,.2f}%",
            },
            "Q30 Bases in Probe Read": {
                "title": "Q30 Bases in Probe Read",
                "description": "Q30 Bases in Probe Read",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:,.2f}%",
            },
            "Estimated UMIs from Genomic DNA": {
                "title": "Estimated UMIs from Genomic DNA",
                "description": "Estimated UMIs from Genomic DNA",
                "scale": "",
                "format": "{:.2f}%",
            },
            "Fraction Reads in Squares Under Tissue": {
                "title": "Fraction Reads in Squares Under Tissue",
                "description": "Fraction Reads in Squares Under Tissue",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:,.2f}%",
            },
            "Mean Genes Under Tissue per Square 2 µm": {
                "title": "Mean Genes Under Tissue per Square 2 µm",
                "description": "Mean Genes Under Tissue per Square 2 µm",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Mean Genes Under Tissue per Bin 8 µm": {
                "title": "Mean Genes Under Tissue per Bin 8 µm",
                "description": "Mean Genes Under Tissue per Bin 8 µm",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Mean Genes Under Tissue per Bin 16 µm": {
                "title": "Mean Genes Under Tissue per Bin 16 µm",
                "description": "Mean Genes Under Tissue per Bin 16 µm",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Median Genes per Cell": {
                "title": "Median Genes per Cell",
                "description": "Median Genes per Cell",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Fraction Reads in Spots Under Tissue": {
                "title": "Fraction Reads in Spots Under Tissue",
                "description": "Fraction Reads in Spots Under Tissue",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:,.2f}%",
            },
            "Number of Spots Under Tissue": {
                "title": "Number of Spots Under Tissue",
                "description": "Number of Spots Under Tissue",
                "scale": "",
                "format": "{:,.0f}",
            },
            "Median Genes per Spot": {
                "title": "Median Genes per Spot",
                "description": "Median Genes per Spot",
                "scale": "",
                "format": "{:,.2f}",
            },
        }
        return table.plot(
            data_by_sample,
            headers,
            pconfig=TableConfig(id="Space Ranger table", title="Space Ranger table: Data Quality"),
        )
