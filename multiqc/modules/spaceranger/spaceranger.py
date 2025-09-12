import logging
from collections import defaultdict
import json
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import TableConfig


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    h/t: https://github.com/MultiQC/MultiQC/blob/main/multiqc/modules/xenium/xenium.py
    Todo
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
                data_by_sample[sample_name] = parsed_data
                self.add_data_source(f, sample_name)

        for f2 in self.find_log_files("spaceranger/count_html", filehandles=True):
            parsed_data_count = self.parse_count_html(f2)
            if parsed_data_count:
                sample_name = parsed_data_count["Sample ID"]
                try:
                    data_by_sample[sample_name] = {**data_by_sample[sample_name], **parsed_data_count}
                except KeyError:
                    data_by_sample[sample_name] = parsed_data_count
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
                    plot=self.spaceranger_transcript_table(data_by_sample))
 
        self.add_section(
                    name="Sequencing Saturation",
                    anchor="Sequencing Saturation",
                    description="Plot of sequencing saturation",
                    helptext="Plot of sequencing saturation",
                    plot=self.add_seq_sat_plot(data_by_sample))
 
        self.add_section(
                    name="Genomic UMIs",
                    anchor="Genomic UMIs",
                    description="Plot of Genomic UMIs",
                    helptext="""
                    Plot of Genomic UMIs
                    """,
                    plot=self.add_gen_umi_plot(data_by_sample))

        self.add_section(name="Genes detected",
                    anchor="Plot",
                    description="Plot of Genes detected in differing sized bins",
                    helptext="""
                    Plot of Genes detected in differing sized bins
                    """,
                    plot=self.add_gene_number_plot(data_by_sample))
     
  
    def add_seq_sat_plot(self, data_by_sample):
        config = {"ylab": "Sequencing Saturation (%)", "cpswitch": False}
        config["id"] = "Space Ranger Sequencing Saturation plot"
        config["title"] = "Space Ranger: Saturation plot"
        
        genes_detected = {
            "Sequencing Saturation": {"color": "#f7a35c", "name": "Sequencing Saturation"}
            }
        return(bargraph.plot(data_by_sample, genes_detected, config))


    def parse_count_html(self, f):
        """
        Space Ranger count report parser
        """

        warnings_data_by_sample: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict(lambda: defaultdict())

        warnings_headers: Dict = dict()
        summary = None

        for line in f["f"]:
            line = line.strip()
            if line.startswith("const data"):
                line = line.replace("const data = ", "")
                summary = json.loads(line)
                if 'summary' in summary.keys():
                    summary = summary["summary"]
                break

        if summary is None:
            logging.error(f"Couldn't find JSON summary data in HTML report, skipping: {f['fn']}")

        sample_name = self.clean_s_name(summary["sample"]["id"], f)

        # List of data collated from different tables in cellranger reports.
        # This is a list of Tuples (metric name, value)

        # For spacreanger 1 & 2:
        if 'summary_tab' in summary.keys():
            software = next(
                iter(x[1] for x in summary["summary_tab"]["pipeline_info_table"]["rows"] if x[0] == "Pipeline Version")
            )

            data_rows = (
                summary['summary_tab']['pipeline_info_table']['rows']
            )

        # For spaceranger 3 & 4:
        if 'tabs' in summary.keys():
            software = next(
                iter(x[1] for x in summary['tabs']['tab_data'][0]['run_summary']['card']['inner']['rows'] if x[0] == "Pipeline Version")
            )

            data_rows = ( 
                    summary['tabs']['tab_data'][0]['run_summary']['card']['inner']['rows']
                )
            
        
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

        return html_dict


    def add_gen_umi_plot(self, data_by_sample):
        config = {"ylab": "Fraction Genomic UMIs (%)", "cpswitch": False}
        
        config["id"] = "Space Ranger Genomic UMIs plot"
        config["title"] = "Space Ranger: Genomic UMIs"
        
        genes_detected = {
            "Estimated UMIs from Genomic DNA": {"color": "#320faf", "name": "Sequencing Saturation"}
            }
        return(bargraph.plot(data_by_sample, genes_detected, config))
    
    
    def add_gene_number_plot(self, data_by_sample):
        config = {"ylab": "Genes detected by bin size"}
        config["id"] = "Space Ranger: Genes detected by bin size"
        config["title"] = "Space Ranger: Genes detected by bin size - beta"
        genes_detected = {
            "Mean Genes Under Tissue per Square 2 µm": {"color": "#20568f", "name": "Mean Genes Under Tissue per Square 2 µm"},
            "Mean Genes Under Tissue per Bin 8 µm": {"color": "#f7a35c", "name": "Mean Genes Under Tissue per Bin 8 µm"},
            "Mean Genes Under Tissue per Bin 16 µm": {"color": "#981919", "name": "Mean Genes Under Tissue per Bin 16 µm"}
            }
        
        return(bargraph.plot(data_by_sample, genes_detected, config))
    
       
    def parse_spaceranger_metrics(self, f):
        """Parse Space Ranger metrics_summary.csv file"""
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

        # Convert numeric fields
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
            "Median Cell Area (μm²)"
            "Median Nucleus Area (μm²)",
            "Maximum Nucleus Diameter (pixels)",
            # visium SD
            'Number of Spots Under Tissue',
            'Number of Reads',
            'Mean Reads per Spot',
            'Mean Reads Under Tissue per Spot',
            'Fraction of Spots Under Tissue',
            'Median Genes per Spot',
            'Median UMI Counts per Spot',
            'Valid Barcodes',
            'Valid UMIs',
            'Sequencing Saturation',
            'Q30 Bases in Barcode',
            'Q30 Bases in RNA Read',
            'Q30 Bases in UMI',
            'Reads Mapped to Genome',
            'Reads Mapped Confidently to Genome',
            'Reads Mapped Confidently to Intergenic Regions',
            'Reads Mapped Confidently to Intronic Regions',
            'Reads Mapped Confidently to Exonic Regions',
            'Reads Mapped Confidently to Transcriptome',
            'Reads Mapped Antisense to Gene',
            'Fraction Reads in Spots Under Tissue',
            'Total Genes Detected',
            # some Visium3
            'Total Genes Detected',
            'Reads Mapped to Genome',
            'Reads Mapped Confidently to Genome',
            'Reads Mapped Confidently to Intergenic Regions',
            'Reads Mapped Confidently to Intronic Regions',
            'Reads Mapped Confidently to Exonic Regions'
        ]

        parsed_metrics = {}
        for field in numeric_fields:
            if field in metrics and metrics[field] != "":
                try:
                    parsed_metrics[field] = float(metrics[field])
                except ValueError:
                    log.warning(f"Could not convert {field}='{metrics[field]}' to float")

        # Keep string fields
        string_fields = [
            "Sample ID"
        ]
        for field in string_fields:
            if field in metrics:
                parsed_metrics[field] = metrics[field]

        # covert some field to percent    
        covert_to_percent = [
            "Sequencing Saturation",
            "Estimated UMIs from Genomic DNA"
        ]
        
        if field in parsed_metrics and parsed_metrics[field] != "":
            for field in covert_to_percent:
                if field in parsed_metrics:
                    parsed_metrics[field] = parsed_metrics[field] * 100
        
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
            "Transcriptome": {
                "title": "Transcriptome",
                "description": "Transcriptome",
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
            },"Sequencing Saturation": {
                "title": "Sequencing Saturation",
                "description": "Sequencing Saturation",
                "scale": "",
                "format": "{:.2f}%",
            },"Reads Mapped to Probe Set": {
                "title": "Reads Mapped to Probe Set",
                "description": "Reads Mapped to Probe Set",
                "modify": lambda x: x * 100.0,
                "scale": "",
                "format": "{:.2f}%",
            },"Reads Mapped Confidently to Probe Set": {
                "title": "Reads Mapped Confidently to Probe Set",
                "description": "Reads Mapped Confidently to Probe Set",
                "modify": lambda x: x * 100.0,
                "scale": "",
                "format": "{:.2f}%",
            },"Reads Mapped Confidently to the Filtered Probe Set": {
                "title": "Reads Mapped Confidently to the Filtered Probe Set",
                "description": "Reads Mapped Confidently to the Filtered Probe Set",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },"Reads Mapped to Genome": {
                "title": "Reads Mapped to Genome",
                "description": "Reads Mapped to Genome",
                "modify": lambda x: x * 100.0,
                "scale": "",
                "format": "{:.2f}%",
            },"Reads Mapped Confidently to Genome": {
                "title": "Reads Mapped Confidently to Genome",
                "description": "Reads Mapped Confidently to Genome",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },"Reads Mapped Confidently to Exonic Regions": {
                "title": "Reads Mapped Confidently to Exonic Regions",
                "description": "Reads Mapped Confidently to Exonic Regions",
                "scale": "",
                "modify": lambda x: x * 100.0,
                "format": "{:.2f}%",
            },"Total Genes Detected": {
                "title": "Genes Detected (Genome)",
                "description": "Genes Detected (Genome)",
                "scale": "",
                "format": "{:,.0f}",
            },"Genes Detected": {
                "title": "Genes Detected (Probe Set)",
                "description": "Genes Detected (Probe Set)",
                "scale": "",
                "format": "{:,.0f}",
            },"Estimated UMIs from Genomic DNA": {
                "title": "Estimated UMIs from Genomic DNA",
                "description": "Estimated UMIs from Genomic DNA",
                "scale": "",
                "format": "{:.2f}%",
            },"Fraction Reads in Squares Under Tissue": {
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
        return table.plot(data_by_sample,headers,pconfig=TableConfig(id="Space Ranger table",title="Space Ranger table: Data Quality"))
