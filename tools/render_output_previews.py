"""Render the first three physical lines of pipeline text outputs as JPEGs.

Uses actual output text, Arial labels and Consolas content. Images and per-trait
.sig filenames are listed without their contents. Requires Pillow and the two font files.
"""

import argparse
import csv
import hashlib
import json
from itertools import islice
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from PIL import Image, ImageDraw, ImageFont


OUTPUTS = {
    1: ['step1.summary.csv', 'step1.PC.csv', 'step1.markers.csv', 'step1.eigenvalues.csv'],
    2: ['step2.all_sig.csv'],
    3: ['step3.all_sig_gene.csv', 'step3.all_unique_snps.csv'],
    4: ['step4.all_sig_l.csv', 'step4.all_locus.csv'],
    5: ['step5.all_locus_cata.csv'],
    6: ['step6.category_summary.csv'],
    7: ['step7.candidate_locus_cata_ALL.csv'],
}
FIGURE_NUMBERS = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 7, 7: 8}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('tests'))
    parser.add_argument('--output', type=Path, default=Path('tests/expected_outputs'))
    parser.add_argument('--font-dir', type=Path, default=Path('C:/Windows/Fonts'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    label = ImageFont.truetype(str(args.font_dir / 'arial.ttf'), 32)
    heading = ImageFont.truetype(str(args.font_dir / 'arial.ttf'), 38)
    mono = ImageFont.truetype(str(args.font_dir / 'consola.ttf'), 22)
    shape_font = ImageFont.truetype(str(args.font_dir / 'consola.ttf'), 18)
    max_width, margin, inset, line_height = 3300, 44, 24, 36
    available = max_width - 2 * (margin + inset)
    measure = ImageDraw.Draw(Image.new('RGB', (1, 1)))
    manifest = {'source_directory': str(args.source.resolve()),
                'fonts': {'filename': 'Arial', 'content': 'Consolas'},
                'tab_width': 4, 'maximum_source_lines': 3,
                'excluded': ['*.sig contents', 'image outputs'], 'steps': []}
    figures = []
    assets = []
    for step, names in OUTPUTS.items():
        panels = []
        for name in names:
            path = args.source / name
            with path.open(encoding='utf-8-sig', newline='') as handle:
                lines = [line.rstrip('\r\n') for line in islice(handle, 3)]
            with path.open(encoding='utf-8-sig', newline='') as handle:
                reader = csv.reader(handle, delimiter='\t')
                header = next(reader)
                column_count = len(header)
                row_count = 0
                for record in reader:
                    if not record:
                        continue
                    if len(record) != column_count:
                        raise ValueError(f'{path}: inconsistent column count')
                    row_count += 1
            shown, clipped = [], []
            for index, line in enumerate(lines):
                expanded = line.expandtabs(4)
                if measure.textlength(expanded, font=mono) > available:
                    end = int(available / measure.textlength('0', font=mono)) - 1
                    expanded = expanded[:end] + '\u2026'
                    assert measure.textlength(expanded, font=mono) <= available
                    clipped.append(index + 1)
                shown.append(expanded)
            panels.append({'filename': name, 'source_lines': lines,
                           'shape': {'rows': row_count, 'columns': column_count},
                           'display_lines': shown, 'clipped_lines': clipped,
                           'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        filename_only = ([p.relative_to(args.source).as_posix()
                          for p in sorted((args.source / 'sig').glob('*.sig'))]
                         if step == 2 else [])
        # Keep native type sizes identical; trim unused horizontal space.
        content_width = max(measure.textlength(line, font=mono)
                            for panel in panels for line in panel['display_lines'])
        width = min(max_width, max(1100, int(content_width) + 2*(margin+inset) + 2))
        # A fixed three-line viewport also makes files with fewer lines recognizable.
        panel_height = 246
        height = 100 + panel_height * len(panels) + 48 * len(filename_only)
        canvas = Image.new('RGB', (width, height), '#141414')
        draw = ImageDraw.Draw(canvas)
        draw.text((margin, 28), f'Step {step}', font=heading, fill='#F5F5F5')
        for index, panel in enumerate(panels):
            y = 100 + panel_height * index
            draw.text((margin, y), panel['filename'], font=label, fill='#F5F5F5')
            top, bottom = y + 48, y + 220
            draw.rectangle((margin, top, width-margin, bottom), fill='#222222',
                           outline='#707070', width=2)
            for row, line in enumerate(panel['display_lines']):
                draw.text((margin+inset, top+17+row*line_height), line, font=mono,
                          fill='#F0F0F0' if row == 0 else '#E6E6E6')
            shape = panel['shape']
            shape_text = f"{shape['rows']:,} rows \u00d7 {shape['columns']} columns"
            draw.text((width-margin-inset, bottom-14), shape_text,
                      font=shape_font, fill='#999999', anchor='rb')
        for index, name in enumerate(filename_only):
            draw.text((margin, 100 + panel_height * len(panels) + 48 * index),
                      name, font=label, fill='#F5F5F5')
        figure_number = FIGURE_NUMBERS[step]
        filename = f'Supplementary Figure S{figure_number}.jpg'
        destination = args.output / filename
        canvas.save(destination, format='JPEG', quality=95, subsampling=0, dpi=(300, 300))
        with Image.open(destination) as check:
            check.verify()
        assets.append(destination)
        manifest['steps'].append({'step': step, 'image': filename,
                                  'figure_number': figure_number,
                                  'title': f'Supplementary Figure S{figure_number}. Example contents of the output files from Step {step}.',
                                  'dimensions': [width, height], 'files': panels,
                                  'filenames_only': filename_only})
        figures.append(f'<figure><h2>Step {step}</h2><a href="{filename}">'
                       f'<img src="{filename}" alt="Step {step} text outputs"></a>'
                       '</figure>')
    manifest_path = args.output / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    gallery = args.output / 'index.html'
    gallery.write_text('''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Pipeline output previews</title><style>
body{margin:32px auto;padding:0 24px;max-width:1600px;background:#141414;color:#eeeeee;font-family:Arial,sans-serif}
figure{margin:40px 0}img{display:block;width:100%;height:auto}a{color:#dddddd}
figcaption,p{color:#b6b6b6;line-height:1.6}h2{font-size:22px}
</style><h1>Pipeline output previews</h1><p>Actual example outputs, first three lines including the header.
Per-phenotype .sig filenames are listed without content; image outputs are excluded. Click an image to inspect it at full resolution.
Tabs use four-character stops; an ellipsis marks right-edge clipping. Full source excerpts are preserved
in <a href="manifest.json">manifest.json</a>.</p>''' + '\n'.join(figures) + '</html>', encoding='utf-8')
    with ZipFile(args.output / 'step_output_screenshots.zip', 'w', ZIP_DEFLATED) as archive:
        for path in assets + [manifest_path, gallery]:
            archive.write(path, path.name)
    print(json.dumps({'images': len(assets), 'files': sum(map(len, OUTPUTS.values())),
                      'output_directory': str(args.output.resolve()),
                      'clipped_files': [p['filename'] for s in manifest['steps']
                                        for p in s['files'] if p['clipped_lines']]}, indent=2))


if __name__ == '__main__':
    main()
