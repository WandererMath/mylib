from Bio.Seq import Seq

def translate(seq):
    dna = Seq(seq)
    x=dna.translate(table=1, to_stop=True)
    return str(x)

import colorsys
def generate_distinct_colors(names):
    """Generate n evenly spaced colors (hex format)."""
    n=len(names)
    colors = []
    for i in range(n):
        hue = i / n  # evenly spaced hues around the color wheel
        lightness = 0.5  # middle brightness
        saturation = 0.65  # medium saturation
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        hex_color = f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}"
        colors.append(hex_color)
    name_to_color = dict(zip(names, colors))
    return name_to_color

def generate_blue_scale(n):
    """Generate n evenly spaced shades of blue (hex format)."""
    colors = []
    for i in range(n):
        # interpolate between light blue and dark blue
        # lightness goes from 0.9 (light) to 0.3 (dark)
        lightness = 0.9 - (0.6 * i / (n - 1)) if n > 1 else 0.5
        hue = 220 / 360  # blue hue
        saturation = 0.8
        import colorsys
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        hex_color = f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}"
        colors.append(hex_color)
    return colors


def generate_red_scale(n):
    """Generate n evenly spaced shades of blue (hex format)."""
    colors = []
    for i in range(n):
        # interpolate between light blue and dark blue
        # lightness goes from 0.9 (light) to 0.3 (dark)
        lightness = 0.9 - (0.6 * i / (n - 1)) if n > 1 else 0.5
        hue = 0 / 360  
        saturation = 0.8
        import colorsys
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        hex_color = f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}"
        colors.append(hex_color)
    return colors



TEMPLATE_1='DATASET_COLORSTRIP\nSEPARATOR TAB\n'

TEMPLATE_2='''TREE_COLORS
SEPARATOR TAB
DATA\n'''


def itol_color_strip(NAME_TO_GENUS:dict, output, GENUS_COLOR_MAPPING=None):
    if GENUS_COLOR_MAPPING is None:
        GENUS_COLOR_MAPPING=generate_distinct_colors(set(NAME_TO_GENUS.values()))
    with open(output, 'w') as f:
        f.write(TEMPLATE_1)
        f.write("DATASET_LABEL")
        for genus in GENUS_COLOR_MAPPING:
            f.write(f"\t{genus}")
        f.write("\n")
        f.write("COLOR")
        for genus in GENUS_COLOR_MAPPING:
            f.write(f"\t{GENUS_COLOR_MAPPING[genus]}")
        f.write("\nDATA\n")

        for protein in NAME_TO_GENUS:
            genus=NAME_TO_GENUS[protein]
            color=GENUS_COLOR_MAPPING[genus]
            f.write(f"{protein}\t{color}\t{genus}\n")  



def itol_color_strip_gradient(SPECIES_PLASMIDS:dict, output):
    with open(output, 'w') as f:
        f.write(TEMPLATE_1)
        f.write("DATASET_LABEL")
        lower=min(SPECIES_PLASMIDS.values())
        red_scale=generate_red_scale(max([v for v in SPECIES_PLASMIDS.values()])+1-lower)
        for i, c in enumerate(red_scale):
            f.write(f"\t{i+lower}")
        f.write("\n")
        f.write("COLOR")
        for c in red_scale:
            f.write(f"\t{c}")
        f.write("\nDATA\n")        
        # breakpoint()
        for species in SPECIES_PLASMIDS:
            num_plasmids=SPECIES_PLASMIDS[species]
            f.write(f"{species}\t{red_scale[num_plasmids-lower]}\t{num_plasmids}\n")
    