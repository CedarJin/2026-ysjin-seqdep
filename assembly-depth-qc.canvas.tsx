import { BarChart, Callout, Card, CardBody, CardHeader, Divider, Grid, H1, H2, LineChart, Stack, Stat, Table, Text } from 'cursor/canvas';

const depthRows = [
  ['10M', '75', '38.30%', '66.84%', '69.14%', '84.99%', '755.86 MB', '893.86 MB', '1.03 GB', '9.57%', '29.98 MB', '54.29 MB', '181.09 MB', '64.49%'],
  ['20M', '75', '38.93%', '70.20%', '71.38%', '87.35%', '1.42 GB', '1.69 GB', '1.97 GB', '9.20%', '46.81 MB', '77.42 MB', '265.16 MB', '64.37%'],
  ['30M', '75', '39.10%', '71.43%', '72.27%', '88.13%', '2.08 GB', '2.47 GB', '2.86 GB', '9.10%', '57.49 MB', '94.07 MB', '327.87 MB', '64.39%'],
  ['40M', '75', '39.19%', '72.06%', '72.74%', '88.77%', '2.71 GB', '3.22 GB', '3.72 GB', '9.18%', '65.44 MB', '106.66 MB', '375.00 MB', '64.49%'],
  ['50M', '75', '39.23%', '72.49%', '73.09%', '89.16%', '3.31 GB', '3.94 GB', '4.56 GB', '9.31%', '72.08 MB', '116.42 MB', '413.38 MB', '64.69%'],
];

const omicRows = [
  ['metaG', '10M', '40', '38.30%', '62.42%', '57.91%', '66.88%'],
  ['metaG', '20M', '40', '38.93%', '65.41%', '60.58%', '70.22%'],
  ['metaG', '30M', '40', '39.10%', '66.27%', '61.67%', '71.45%'],
  ['metaG', '40M', '40', '39.19%', '66.74%', '62.28%', '72.08%'],
  ['metaG', '50M', '40', '39.23%', '67.04%', '62.65%', '72.55%'],
  ['metaT', '10M', '35', '77.04%', '82.15%', '81.97%', '84.99%'],
  ['metaT', '20M', '35', '78.66%', '84.03%', '83.72%', '87.35%'],
  ['metaT', '30M', '35', '77.97%', '84.60%', '84.37%', '88.13%'],
  ['metaT', '40M', '35', '79.90%', '85.01%', '84.69%', '88.77%'],
  ['metaT', '50M', '35', '80.18%', '85.19%', '85.01%', '89.16%'],
];

const bamVariationRows = [
  ['metaT', 'MT0001', '10M', '0.38%', '755.86 MB', '762.78 MB', '77.04-78.76%'],
  ['metaT', 'MT0001', '50M', '0.25%', '3.31 GB', '3.33 GB', '80.18-81.23%'],
  ['metaT', 'MT0004', '40M', '0.25%', '2.93 GB', '2.95 GB', '84.73-85.30%'],
  ['metaT', 'MT0004', '20M', '0.24%', '1.56 GB', '1.57 GB', '83.51-84.36%'],
  ['metaT', 'MT0007', '50M', '0.23%', '3.83 GB', '3.86 GB', '83.41-83.96%'],
];

const contigVariationRows = [
  ['metaG', 'MG0005', '10M', '0.71%', '178.18 MB', '181.09 MB', '62.56-62.88%'],
  ['metaG', 'MG0004', '50M', '0.53%', '287.49 MB', '291.13 MB', '62.16-62.23%'],
  ['metaT', 'MT0005', '20M', '0.26%', '108.04 MB', '108.70 MB', '86.97-87.35%'],
  ['metaT', 'MT0006', '50M', '0.23%', '87.24 MB', '87.73 MB', '82.83-83.56%'],
  ['metaT', 'MT0004', '10M', '0.21%', '49.34 MB', '49.61 MB', '81.27-82.41%'],
];

const lowAlignmentRows = [
  ['10M', 'metaG/MG0007', '38.30-38.51%', 'all five seeds'],
  ['20M', 'metaG/MG0007', '38.93-39.00%', 'all five seeds'],
  ['30M', 'metaG/MG0007', '39.10-39.23%', 'all five seeds'],
  ['40M', 'metaG/MG0007', '39.19-39.37%', 'all five seeds'],
  ['50M', 'metaG/MG0007', '39.23-39.35%', 'all five seeds'],
];

const depths = ['10M', '20M', '30M', '40M', '50M'];

export default function AssemblyDepthQc() {
  return (
    <Stack gap={18}>
      <H1>Assembly / Alignment Depth QC</H1>
      <Text tone="secondary">
        Checked 375 final alignments grouped by sub sequencing depth, sample, omic, and seed.
      </Text>

      <Grid columns={4} gap={14}>
        <Stat value="375" label="BAM + BAI outputs checked" tone="success" />
        <Stat value="0" label="missing or empty core outputs" tone="success" />
        <Stat value="0.38%" label="largest BAM seed CV" tone="success" />
        <Stat value="0.71%" label="largest final-contig seed CV" tone="success" />
      </Grid>

      <Callout tone="success" title="Size consistency looks good">
        For every omic/sample/depth group, all five seeds are present. Seed-to-seed size variation is very small:
        BAM CV max is 0.38%, and final-contig CV max is 0.71%.
      </Callout>

      <Grid columns="1.2fr 1fr" gap={16}>
        <Card>
          <CardHeader>Median alignment rate by depth</CardHeader>
          <CardBody>
            <LineChart
              categories={depths}
              valueSuffix="%"
              height={220}
              series={[
                { name: 'All outputs median', data: [66.84, 70.20, 71.43, 72.06, 72.49], tone: 'info' },
                { name: 'metaG median', data: [62.42, 65.41, 66.27, 66.74, 67.04] },
                { name: 'metaT median', data: [82.15, 84.03, 84.60, 85.01, 85.19], tone: 'success' },
              ]}
            />
          </CardBody>
        </Card>

        <Card>
          <CardHeader>Median BAM size by depth</CardHeader>
          <CardBody>
            <BarChart
              categories={depths}
              valueSuffix=" GB"
              height={220}
              series={[{ name: 'Median BAM size', data: [0.87, 1.69, 2.47, 3.22, 3.94], tone: 'info' }]}
            />
          </CardBody>
        </Card>
      </Grid>

      <Divider />
      <H2>Depth Summary</H2>
      <Table
        headers={['Depth', 'N', 'Rate min', 'Rate median', 'Rate mean', 'Rate max', 'BAM min', 'BAM median', 'BAM max', 'BAM CV', 'Contig min', 'Contig median', 'Contig max', 'Contig CV']}
        rows={depthRows}
        columnAlign={['left', 'right', 'right', 'right', 'right', 'right', 'right', 'right', 'right', 'right', 'right', 'right', 'right', 'right']}
        striped
      />

      <H2>Alignment Rate by Omic and Depth</H2>
      <Table
        headers={['Omic', 'Depth', 'N', 'Min', 'Median', 'Mean', 'Max']}
        rows={omicRows}
        columnAlign={['left', 'left', 'right', 'right', 'right', 'right', 'right']}
        striped
      />

      <Grid columns={2} gap={16}>
        <Card>
          <CardHeader>Largest BAM size variation across seeds</CardHeader>
          <CardBody>
            <Table
              headers={['Omic', 'Sample', 'Depth', 'CV', 'Min', 'Max', 'Rate range']}
              rows={bamVariationRows}
              columnAlign={['left', 'left', 'left', 'right', 'right', 'right', 'right']}
              framed={false}
              striped
            />
          </CardBody>
        </Card>

        <Card>
          <CardHeader>Largest final-contig size variation across seeds</CardHeader>
          <CardBody>
            <Table
              headers={['Omic', 'Sample', 'Depth', 'CV', 'Min', 'Max', 'Rate range']}
              rows={contigVariationRows}
              columnAlign={['left', 'left', 'left', 'right', 'right', 'right', 'right']}
              framed={false}
              striped
            />
          </CardBody>
        </Card>
      </Grid>

      <H2>Lowest Alignment Rates</H2>
      <Text tone="secondary">
        The consistently lowest alignment group is metaG/MG0007 at every depth. This looks biological or sample-specific rather than a failed job, because all five seeds are present and close to each other.
      </Text>
      <Table
        headers={['Depth', 'Group', 'Alignment rate range', 'Seeds affected']}
        rows={lowAlignmentRows}
        columnAlign={['left', 'left', 'right', 'left']}
        rowTone={['warning', 'warning', 'warning', 'warning', 'warning']}
      />
    </Stack>
  );
}
