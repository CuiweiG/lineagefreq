WHO FluNet data requires interactive download from:
https://www.who.int/tools/flunet

The bulk download URLs (frontdoor-*.azurefd.net and apps.who.int/flumart)
either return 404 or redirect to an interactive HTML page.

To obtain this data manually:
1. Visit https://www.who.int/tools/flunet
2. Select parameters (region, date range, virus type)
3. Click "Download CSV"
4. Save the file as who_flunet.csv in this directory

The lineagefreq validation uses the included simulated influenza_h3n2
dataset as an alternative for multi-pathogen demonstration.
