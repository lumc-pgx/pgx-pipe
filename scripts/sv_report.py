from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot 
from bokeh.models import HoverTool, BoxAnnotation
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.models.widgets import Div

import gffutils

import gene_viz
from gene_viz.gene_viz import zooming_ticker
from gene_viz.features import Transcript, Exon, CDS
from gene_viz.utils import transcripts_from_gffutils


class segment(object):
    """
    A representation of an aligned sequence fragment
    """
    def __init__(self, ident, start, end, strand):
        self.id = ident
        self.start = start
        self.end = end
        self.strand = strand
    
    def __str__(self):
        return "< " + "\t".join([self.id, self.start, self.end, self.strand]) + " >"
    
    def __repr__(self):
        return self.__str__()

        
class alignment(object):
    """
    A representation of a contigious alignment of two sequences
    """
    def __init__(self, query, hit):
        self.query = query
        self.hit = hit
    
    def __str__(self):
        return "< " + "\t".join([str(self.query), str(self.hit)]) + " >"
    
    def __repr__(self):
        return self.__str__()
        
    @classmethod
    def from_string(cls, string):
        """
        Construct an alignment object from a string containing a textual representation
        """
        fields = string.strip().split('\t')
        assert len(fields) == 11
        query = segment(fields[0], int(fields[2]), int(fields[3]), fields[4])
        hit = segment(fields[5], int(fields[7]), int(fields[8]), fields[9])
        if hit.strand == "-":
            hit.start, hit.end = hit.end, hit.start
        return cls(query, hit)
        
        
def load_alignments(path):
    """
    Read alignments from the file specified by 'path'
    """
    with open(path, "r") as infile:
        for line in infile:
            if line.startswith('#') or line.strip() == "":
                continue
            yield alignment.from_string(line)


# tooltip definition
hover = HoverTool(tooltips=[
    ("Query", " @y0-@y1"),
    ("Reference", " @x0-@x1"),
])


def plot_sv_alleles(annotations, aligned_regions_file):
    alignments = list(load_alignments(aligned_regions_file))
    query_ids = sorted(list(set(a.query.id for a in alignments)))

    plots = []
    for query_id in query_ids:
        alns = [a for a in alignments if a.query.id == query_id]
        x0 = [a.hit.start for a in alns]
        x1 = [a.hit.end for a in alns]
        y0 = [a.query.start for a in alns]
        y1 = [a.query.end for a in alns]
    
        min_x = min(x0 + x1)
        max_x = max(x0 + x1)
        x_range = (min_x, max_x)
        width = x_range[1] - x_range[0]
        pad = 0.1 * width
        x_range = (int(min_x - pad), int(max_x + pad))
        colors = ["red" if a.hit.strand == '-' else "blue" for a in alns]
        fig = figure(width=800, height=300, x_range=x_range)
        fig.segment(x0, y0, x1, y1, line_color=colors, line_width=5, 
                    line_alpha=0.75, line_cap="butt")
        fig.title.text = query_id
        fig.yaxis.axis_label = "Amplicon"
        fig.xaxis.visible = False
        fig.yaxis.ticker = zooming_ticker()
        fig.xaxis.formatter = NumeralTickFormatter(format="0,0")
        fig.add_tools(hover)
        plots.append(fig)
    
        transcripts = transcripts_from_gffutils(annotations, 
            next(a.hit.id for a in alns), x_range[0], x_range[1])
            
        fig2 = gene_viz.GenePlot(
            dict(
                row_height=15,
                exon_color="Black",
                pack=True,
                label_horiz_position="center",
                label_vert_position="above",
                label_justify="center",
                label_offset=(0, -4),
                label_font_size="6pt",
                intron_width_percent=0.001
            )
        )
        
        fig2.link_x_range(fig)
        fig2.figure.xaxis.axis_label = next(a.hit.id for a in alns)
        fig2.figure.yaxis.visible = False
        fig2.figure.background_fill_color = "black"
        fig2.figure.background_fill_alpha = 0.15
        fig2.transcripts = transcripts
        box = BoxAnnotation(left=min_x, right=max_x, fill_color='white', fill_alpha=1, level="underlay")
        fig2.figure.add_layout(box)
        
        for a in alns:
            box = BoxAnnotation(left=min([a.hit.start, a.hit.end]), right=max([a.hit.start, a.hit.end]),
                                fill_color="red" if a.hit.strand == '-' else "blue", fill_alpha=0.25, level="underlay")
            fig2.figure.add_layout(box)
            
        fig2.update()
        plots.append(fig2.figure)
    
    if len(plots) > 0:
        save(gridplot(plots, ncols=1, toolbar_location="right", toolbar_options=dict(logo=None)))
    else:
        save(Div(text="No Data"))


# set bokeh to write to file
output_file(snakemake.output[0], title=snakemake.wildcards.barcode)

# connect to the annotation database
db = gffutils.FeatureDB(snakemake.input.annotations)

# make the figures
plot_sv_alleles(db, snakemake.input.alignments)
