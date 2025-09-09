import logging
import pandas as pd

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, box, linegraph, scatter, table
from multiqc.plots.table_object import ColumnDict, TableConfig
from multiqc.utils import mqc_colour


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
   
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Space Ranger",
            anchor="spaceranger",
            href=["https://www.10xgenomics.com/support/software/space-ranger/latest"],
            info="Space Ranger is a set of analysis pipelines that process 10x Genomics Visium data with brightfield or fluorescence microscope images, allowing users to map the whole transcriptome in a variety of tissues",
        )

        data_by_sample = {}
        for f in self.find_log_files("spaceranger", filehandles=True):
            parsed_data = self.parse_spaceranger_metrics(f)
            if parsed_data:
                sample_name = parsed_data["Sample ID"]
                data_by_sample[sample_name] = parsed_data
                self.add_data_source(f, sample_name)
        
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
                    helptext="""
                    Plot of sequencing saturation
                    """,
                    plot=self.add_set_sat_plot(data_by_sample))
 
        self.add_section(
                    name="Genes detected",
                    anchor="Plot",
                    description="Plot of Genes detected in differing sized bins",
                    helptext="""
                    Plot of Genes detected in differing sized bins
                    """,
                    plot=self.add_gene_number_plot(data_by_sample))
     
     
    def add_set_sat_plot(self, data_by_sample):
        config = {"ylab": "Sequencing Saturation (%)"}
        
        config["id"] = "Space Ranger Sequencing Saturation plot"
        config["title"] = "Space Ranger: Saturation plot"
        
        genes_detected = {
            "Sequencing Saturation": {"color": "#f7a35c", "name": "Sequencing Saturation"}
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
            'Total Genes Detected'
            # some Visium3'
            'Reads Mapped to Genome',
            'Reads Mapped Confidently to Genome',
            'Reads Mapped Confidently to Intergenic Regions',
            'Reads Mapped Confidently to Intronic Regions',
            'Reads Mapped Confidently to Exonic Regions'
        ]
        # do duplicates even matter? idk
        #numeric_fields = list(set(numeric_fields))

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
                
        
        parsed_metrics["Sequencing Saturation"] = parsed_metrics["Sequencing Saturation"] * 100
        return parsed_metrics

    def spaceranger_general_stats_table(self, data_by_sample):
        """Add key Xenium metrics to the general statistics table"""
        headers: Dict[str, Dict[str, Any]] = {
            "Number of Reads": {
                "title": "Number of Reads",
                "description": "Total number of reads sequenced",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
             "Reads Mapped to Probe Set": {
                "title": "Reads Mapped to Probe Set",
                "description": "Reads Mapped to Probe Set",
                "scale": "Purples",
                "format": "{:,.2f}",
            },
             "Reads Mapped to Genome": {
                "title": "Reads Mapped to Genome",
                "description": "Reads Mapped to Genome",
                "scale": "Purples",
                "format": "{:,.2f}",
             },"Genes Detected": {
                "title": "Genes Detected",
                "description": "Genes Detected",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
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
                "scale": "",
                "format": "{:,.2f}",
            },"Reads Mapped Confidently to Probe Set": {
                "title": "Reads Mapped Confidently to Probe Set",
                "description": "Reads Mapped Confidently to Probe Set",
                "scale": "",
                "format": "{:,.2f}",
            },"Reads Mapped Confidently to the Filtered Probe Set": {
                "title": "Reads Mapped Confidently to the Filtered Probe Set",
                "description": "Reads Mapped Confidently to the Filtered Probe Set",
                "scale": "",
                "format": "{:,.2f}",
            },"Reads Mapped to Genome": {
                "title": "Reads Mapped to Genome",
                "description": "Reads Mapped to Genome",
                "scale": "",
                "format": "{:,.2f}",
            },"Reads Mapped Confidently to Genome": {
                "title": "Reads Mapped Confidently to Genome",
                "description": "Reads Mapped Confidently to Genome",
                "scale": "",
                "format": "{:,.2f}",
            },"Reads Mapped Confidently to Exonic Regions": {
                "title": "Reads Mapped Confidently to Exonic Regions",
                "description": "Reads Mapped Confidently to Exonic Regions",
                "scale": "",
                "format": "{:,.2f}",
            },"Genes Detected": {
                "title": "Genes Detected",
                "description": "Genes Detected",
                "scale": "",
                "format": "{:,.0f}",
            },"Estimated UMIs from Genomic DNA": {
                "title": "Estimated UMIs from Genomic DNA",
                "description": "Estimated UMIs from Genomic DNA",
                "scale": "",
                "format": "{:,.0f}",
            },"Fraction Reads in Squares Under Tissue": {
                "title": "Fraction Reads in Squares Under Tissue",
                "description": "Fraction Reads in Squares Under Tissue",
                "scale": "",
                "format": "{:,.2f}",
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
