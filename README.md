# Final-Project
# Basic Description
My final project is about breast cancer in Korea, and I collect the data from the website https://dcc.icgc.org/releases/current/Projects/BRCA-KR and gather the data I want.
In my data, all the samples are from young women and they are donors. I want to compare the survival state of the patients whose relapsed type are distant recurrence/metastasis, which can present the survival rate of this relapsed type patients. Thus, in my file, you can see the these data there. 
# Milestone 1 11.13
I'm not sure my data is enough, and in the next seven days, I will move foward to collect data and try to study the raw data of patients' gene if they have some impact in patients' relapsed. Also, I'm thinking about the gene expression also can impact the relapsed rate. Thus, I'm also thinking about use which function and method could to reflect this issue clearly. I will try to use heat map that can reflect the replased rate in the data.
# Milestone 2 11.19
I have gotten my data in R. Next, I will use bar graph to see the survival rate of BRCA and also creat the shiny.io app that allow person to choose the coloring by points of my data. I also conduct pca but I need to think about conduct pca for what. Also, there is a problem that the website of bioinform.io cannot be opened so I cannot see the contents that we learned. 
# Milestone 3 11.27
shiny.io app and pca and all the data anylsis could be done. And figure out the issue in my anylsis. To prepare final presentation. 
# issue
I need to extract the sample_id, gene_id and total_raw_data into a new dataframe. and the sample_id will be the header and gene_id will be the very left colnume and use normalized_read_count to fill the dataframe. After I research the results of the website, I still don't know how to use sample_id to be the header. Now I just know to extract the 3 colnume from the file. And if I use a gene_id to do the sort in the execl, how can I sort other gene_id. However, I will try to use python to figure out this issue if I can do that successfully.
