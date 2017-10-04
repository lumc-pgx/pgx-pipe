import os
from collections import OrderedDict
import datetime
import yaml
from yaml.representer import Representer
from jinja2 import Environment, FileSystemLoader

j2_env = Environment(loader=FileSystemLoader(os.path.dirname(snakemake.input.template)),
                     trim_blocks=True)

j2_template = j2_env.get_template(os.path.basename(snakemake.input.template))

with open(snakemake.output[0], "w") as outfile:
    yaml.add_representer(OrderedDict, Representer.represent_dict)
    print(
        j2_template.render(
            title="Pharmacogenomics Pipeline Configuration Report",
            config = yaml.dump(snakemake.config, default_flow_style=False),
            timestamp="{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now())
        ),
        file=outfile
    )
