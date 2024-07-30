---
title: "ELISA-MS Tool Documentation"
output: md_document #Type of document for output
knit: (knitr::knit)
---




### <b>About ELISA-MS Tool</b><a name = "top"></a>

<p style="text-align: justify;">

The ELISA-MS tool is intended to perform a predictive coverage analysis of new HCP MS analysis to assist the project's principal investigator in making an informed decision to proceed with their ELISA-MS analysis using a specific kit.

</p>


#### Quick-access

1.  [File input](#A)
    -   [Protein Acession Numbers File Type](#A_a)
    -   [Organism Selection](#A_b)
2.  [Analysis](#B)
    -   [Coverage Table](#B_a)
    -   [Coverage Plot](#B_b)
    -   [Venn Diagrams](#B_c)
    -   [PI vs Mw](#B_d)
3.  [Special Features](#C)
    -   [Filtration](#C_a)
    -   [Download](#C_b)

------------------------------------------------------------------------

### <b>1. File input</b><a name = "A"></a>


<p style="text-align: justify;">

The file input section allows users to upload and specify the type of input data they wish to use in the analysis. Users
can choose between different input types such as FASTA files, Excel files, or a list of accession numbers. Each input type has its own specific requirements and functionalities.

</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b>Protein Acession Numbers File Type</b><a name = "A_a"></a>

<p style="text-align: justify;">

<b>Input Type Selection</b>
Users can select the type of input data using the "Input Type" Buttons:
<ol>
<li>FASTA File Input</li>
 - When "FASTA File" is selected, a file input widget appears, allowing users to upload a .fasta file. The file must be in FASTA format, which is a standard text-based format for representing peptide sequences. Currently specific terminologies have to present at the FASTA files at the names of each peptide sequences, specifically: CRIGR, ECOBD, SPOFR
<li>Excel File Input</li>
 - If "Excel File" is selected, the ELISA-MS Tool is expecting a standard Alphalyse SA file. The tool will retrieve the data frame from the *"HCP"* sheet. It is necessary that this sheet has *Protein* and *ppm* column, since some of the functionalities rely on these columns
<li>Accession Number Input</li>
 - Selecting "List of Accession Number" will provoke to a text area to appear for the user to enter accession numbers.
 The user can input multiple accession numbers, separated by commas or newlines (copy pasted from an excel file for example)
</ol>
</p>

**Warning:**
<p style="text-align: justify;">
<!--html_preserve--><i class="fas fa-triangle-exclamation fa-solid" role="presentation" aria-label="triangle-exclamation icon" style="color: #ffc107;"></i><!--/html_preserve--> 
Please ensure that all FASTA file headers include the correct terminologies (CRIGR, ECOBD, SPOFR). Missing or incorrect headers may result in errors or incomplete data processing. For Excel files, verify that the "HCP" sheet contains both *Protein* and *ppm* columns to ensure accurate analysis. If a list of accession numbers separated by commas is provided as input please ensure that all proteins accession numbers are separated by commas.
</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b>Organism Selection</b><a name = "A_b"></a>

<p style="text-align: justify;">

In addition to selecting the input type, users can specify the organism related to their data using the "Organism". This feature was added to filtered out the data of other organism which might not be of interest in the study at hand. Since the ELISA-MS tool is supported by ELISA-MS database, currently we have a limitation of choosing 3 organisms to perform a predictive analysis: *Escherichia coli*, CHO, Spodoptera frugiperda

</p>

<p style="text-align: justify;">
<!--html_preserve--><i class="fas fa-triangle-exclamation fa-solid" role="presentation" aria-label="triangle-exclamation icon" style="color: #ffc107;"></i><!--/html_preserve--> 
By default, *E. coli* is select as the filtering organism. Thus, is essential to change it according to the data at hand. The use of data belonging to other expression models like HEK most likely won't give accurate results. 
</p>


[Move to the top](#top)


------------------------------------------------------------------------

### <b>2. Analysis</b><a name = "B"></a>

<p style="text-align: justify;">

The Analysis section of the tool provides users with a comprehensive set of visual and data analysis options to explore and interpret the input data. After submitting UniProt accession numbers and selecting the relevant organism, the tool processes the data to generate several analytical outputs, each accessible via different tabs in the main panel.

</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b>Coverage Table</b><a name = "B_a"></a>

<p style="text-align: justify;">

The Coverage Table tab presents a detailed view of the protein data, including coverage information for each protein. This table can be searched interactively, allowing users to focus on specific proteins of interest. The table includes columns for protein identifiers, such as UniProt accession numbers.
- <b>Exporting Data:</b> There is an option to download the coverage table as a cstom Excel file, which includes additional metadata such as PI and Mw, as well as kit individual coverage (raw number and percentage)

</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b> Coverage Plot</b><a name = "B_b"></a>

<p style="text-align: justify;">

The Coverage Plot tab features an Upset Plot visualization, created using the 'UpSetR' package, which provides an insightful overview of the intersections between different data sets (e.g., different kits). This plot helps users identify the number unique and shared proteins across various combination of kits.

</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b>Venn Diagrams</b><a name = "B_c"></a>

<p style="text-align: justify;">

In the Venn Diagrams tab, users can visualize the overlap between up to five different datasets. This visualization is particularly useful for understanding the relationships and shared proteins across multiple kits.
- <b>Comparison Selection:</b> Users can select which kits to compare, allowing for tailored analysis based on their specific research questions


</p>

**Warning:**
<p style="text-align: justify;">
<!--html_preserve--><i class="fas fa-triangle-exclamation fa-solid" role="presentation" aria-label="triangle-exclamation icon" style="color: #ffc107;"></i><!--/html_preserve--> 
The user must have at least 2 kits selected in their options for a Venn Diagram to be displayed. Also, a maximum of 5 options can be selected. Otherwise, the Venn Diagram will become overly complex and extremely hard to understand. if more than 5 kits are selected no Diagram will be displayed on the tab
</p>

[Move to the top](#top)

------------------------------------------------------------------------
#### <b>PI vs Mw</b><a name = "B_d"></a>

<p style="text-align: justify;">

The PI vs. MW Plot tab allows users to explore the distribution of proteins based on their isoelectric point (pI) and molecular weight (MW).

- <b>Data Selection:</b> Users can choose up to two kits for comparison, helping to identify proteins unique to each kit or shared between them based on their PI and Mw


</p>

[Move to the top](#top)

------------------------------------------------------------------------
### <b>3. Special Features</b><a name = "C"></a>

<p style="text-align: justify;">

The Special Features section highlights the advanced functionalities of the tool, enabling users to fine-tune their data analysis and customize the output according to their specific needs.

</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b>Filtration</b><a name = "C_a"></a>

<p style="text-align: justify;">

The filtration feature enables users to refine the data displayed in all ELISA-MS tools, from the Coverage Table to the PI vs Mw plot, based on their HCPs' expression levels (parts per million (ppm) expression values). This functionality is critical for identifying the most relevant proteins in a sample.

- <b>Raw Number of HCPs:</b> Users can filter the data to show only the highest-expressing HCPs in a sample. For example, if you set a filter to show the top 100 HCPs, the tool will only show the ones with the highest expression levels. This helps to narrow down the list to the most important proteins.

- <b>Percentage-Based Filtering:</b>This option allows users to filter the data based on a percentage threshold, such
as the top 25% of expressing HCPs. This is useful for comparative studies, ensuring that only the most highly expressed proteins are considered.

- <b>ppm-Based Filtering:</b> For more specific needs, users can filter proteins based solely on their ppm expression values.

These filtration settings significantly impact the data visualizations and the content of the downloadable files, allowing users to focus on the most pertinent data.
</p>

[Move to the top](#top)

------------------------------------------------------------------------

#### <b>Download</b><a name = "C_b"></a>

<p style="text-align: justify;">

The download feature offers extensive customization options, allowing users to save the generated visualizations and data tables in various formats and sizes. This flexibility is vital for documentation, presentation, and further analysis.

- <b> Customizable Dimensions:</b> Users can specify the height and width of the visualizations before downloading them. This feature ensures that the graphics are suitable for different uses, such as presentations, publications, or detailed reports.
-<b> Multiple File Formats:</b> The tool supports downloading the visualizations in multiple formats, including PDF, PNG, and SVG. Each format serves different purposes:
 - <b>PDF:</b> Ideal for high-quality prints and official documentation.
 - <b>PNG:</b> Suitable for web use and presentations due to its wide compatibility.
 - <b>SVG:</b> A scalable vector format that maintains quality at any size, perfect for graphic design and detailed editing.
</p>

[Move to the top](#top)

------------------------------------------------------------------------
