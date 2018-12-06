# Final-Project
# Basic Description
My final project is about breast cancer in Korea, and I collect the data from the website https://dcc.icgc.org/releases/current/Projects/BRCA-KR and gather the data I want.
In my data, all the samples are from young women and they are donors. I want to compare the survival state of the patients whose relapsed type are distant recurrence/metastasis, which can present the survival rate of this relapsed type patients. Thus, in my file, you can see the these data there. 

Milestone 1 11.13
To move foward to collect data and try to study the raw data of patients' gene if they have some impact in patients' relapsed. Also, I'm thinking about the gene expression also can impact the relapsed rate. Thus, I'm also thinking about use which function and method could to reflect this issue clearly. I will try to use heat map that can reflect the replased rate in the data.

Milestone 2 11.19
I have gotten my data in R. Next, I will use bar graph to see the survival rate of BRCA and also creat the shiny.io app that allow person to choose the coloring by points of my data. I also conduct pca but I need to think about conduct pca for what. 

Milestone 3 11.27
shiny.io app and pca and all the data anylsis could be done. And figure out the issue in my anylsis. To prepare final presentation. 

# The Major Milestone
To creat pca to analysis all the data.

# The First Milestone
To collect the data and put the data in R and creat the dataframe.
```{r}
library(ggplot2)
library(dplyr)
library(plotly)
setwd(dir = '~/Desktop/510/final_project')
samples <- read.csv('samples_infor.csv',header = TRUE, sep = ",", 
                          quote = "\"", dec = ".", fill = TRUE)
gene <- read.csv('gene_.csv',header = TRUE, sep = ",", 
                          quote = "\"", dec = ".", fill = TRUE)
```
The issue is my data is not in format that let me use to analysis.

# The Second Milestone
I use python to make my data into a right format. icgc_Sample_id is my header and gene_id is very left column and the value is normalization_read_count.

```{r}
import pandas as pd
import numpy as np

if __name__ == '__main__':
	data = pd.read_csv('exp_seq.BRCA-KR.tsv', sep='\t')

	n_samples = data.shape[0]
	sample_ids = data['icgc_sample_id'].unique()
	gene_ids = data['gene_id'].unique()
	col_names = np.append('gene_id', sample_ids)
	cols = col_names.shape[0]
	rows = gene_ids.shape[0]

	# np.zeros((rows, cols))
	new_data = pd.DataFrame(np.zeros((rows, cols)), columns=col_names, index=gene_ids)
	new_data = new_data.astype('str')
	# print(new_data)

	for r in range(n_samples):
	    if r % 1000 == 0:
	        print(r)
	    col = data.iloc[r]['icgc_sample_id']
	    gene_id = data.iloc[r]['gene_id']
	#     print(gene_id)
	    new_data.loc[gene_id][col] = data.iloc[r]['normalized_read_count']
	    new_data.loc[gene_id]['gene_id'] = gene_id

	new_data.to_csv('data.csv', sep='\t', index=False)
```
Rename the data.csv to exp_read_count.csv
Next, I put data in R into dataframe, and my data need use seq='\t' instead of seq=','
```{r} 
library(ggplot2)
library(reshape2)
library('RColorBrewer')
library(dplyr)
library('dendextend')
library(plotly)

#put data into dataframe
setwd(dir = '~/Desktop/510/final_project') #open the directory
samples <- read.csv('samples_infor.csv',header = TRUE, sep = ",",quote = "\"", dec = ".", fill = TRUE) #read the csv file in R
gene_data <- read.csv('exp_read_count.csv',header = TRUE, sep = "\t",quote = "\"", dec = ".", fill = TRUE) #this file need use "\t" to seperate
```
Next I analysis the donors' vital status and the disease_status_last_followup to make two histogram. And I obtain the data that BRCA is a high complete remission rate and suvival rate disease.
```{r}
#analysis the samples data into the histogram
ggplot(data = samples) +  geom_bar(mapping = aes(x = donor_vital_status), width = 0.5)
ggplot(data = samples) +  geom_bar(mapping = aes(x = disease_status_last_followup, fill = donor_vital_status), width = 0.5)
```
Next I need to creat a heatmap to compare 50 samples with their gene_id.

# The Third Milestone
I created the heatmap using my gene data. However, the issue is my heatmap seems like just has one color which I think is wrong.
```{r}
#there are 50 samples in the dataset and compare them with the gene_id
melted_gene<- melt(gene_data) # the dataset will become two variable column and one value column
ggplot(data = melted_gene, aes(x=melted_gene$gene_id, y=melted_gene$variable, fill=value)) + geom_tile() #the heatmap result
```
Then I create a 50 by 50 matrix of correlation values. A function to help me with this would be “cor” in R, and I can apply it to genes. Also I reshape the data as a list of pairs to make easier to compute. This is core to the concept of “melt”, which turns a square matrix, into minimal (in this case pairwise) compoents. Finally, I put this melted_gene into a heatmap to get a better view.

```{r}
#there are 50 samples in the dataset and compare them
corr_gene <- cor(gene_data[,2:51])
melted_gene<- melt(corr_gene) # the dataset will become two variable column and one value column
head(melted_gene)
p<-ggplot(melted_gene , aes(x = Var1, y = Var2)) + geom_raster(aes(fill = value)) + scale_fill_gradient2(low="navy", mid="white", high="red", midpoint=0.5) + theme( plot.title = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank())
ggplotly(p)
#make a linkage map, and this can be plotted as a dendrogram
clusters <- hclust(dist(gene_data_analysis[,10:11]))
dend <- as.dendrogram(clusters)
dend <- color_branches(dend, k=3)
par(cex=0.5) # reduces font
plot(dend)
```
Next I will creat pca and shiny.io app

# The final major
I creat pca below:

```{r}
pca <- prcomp(gene_data[,2:51],scale = TRUE)
pcadf<- data.frame(pca$rotation)
gene_data_analysis <- cbind(samples,pcadf)
write.csv(gene_data_analysis,"gene_sample_pca.csv")
plot_ly(data = gene_data_analysis, x = ~PC1, y = ~PC2, color = gene_data_analysis$icgc_donor_id )
plot_ly(data = pcadf, x = ~PC1, y = ~PC2, text = rownames(pcadf))
plot_ly(pcadf, x = ~PC1, y = ~PC2, z = ~PC3, color = ~PC4, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
 layout(scene = list(xaxis = list(title = 'PC1'),
                     yaxis = list(title = 'PC2'),
                     zaxis = list(title = 'PC3')))
```

My final project file-finalproject.Rmd has been already uploaded into website:
http://rpubs.com/vidia/FinalProjectBRCA

Also,I creat a data explorer using my gene_sample_pca file. And the publish website as below:
https://xianghuildataexplorer.shinyapps.io/samples/

Finally, I also put my pca into shinyapp.io:https://xianghuildataexplorer.shinyapps.io/PCAshinyapp/
However, there is an error occured when I published. Thus, the code is as below:

```{r}
library(shiny)
library(ggplot2)
library(RColorBrewer)
library(plotly)

dataset <- read.csv('gene_sample_pca.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
headerNames=colnames(dataset)

ui <- fluidPage(
  pageWithSidebar(
    
    headerPanel("PCA Explorer"),
    
    sidebarPanel(
      selectInput('x', 'X', c("None"=FALSE,headerNames),headerNames[11]),
      selectInput('y', 'Y', c("None"=FALSE,headerNames),headerNames[12]),
      selectInput('z', 'z', c("None"=FALSE,headerNames),headerNames[13])
    ),
    
    mainPanel(
      plotlyOutput('plot')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot <- renderPlotly({
    plot_ly(dataset, x = ~get(input$x), y = ~get(input$y), z = ~get(input$z), color = dataset$icgc_donor_id, colors = c('#BF382A', '#0C4B8E')) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')))
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
```
I can run app successful but uploading into the website

# The Completed Final Project Code

## R markdown code part

Code of final project

```{r} 
library(ggplot2)
library(reshape2)
library('RColorBrewer')
library(dplyr)
library('dendextend')
library(plotly)

#put data into dataframe
setwd(dir = '~/Desktop/510/final_project') #open the directory
samples <- read.csv('samples_infor.csv',header = TRUE, sep = ",",quote = "\"", dec = ".", fill = TRUE) #read the csv file in R
gene_data <- read.csv('exp_read_count.csv',header = TRUE, sep = "\t",quote = "\"", dec = ".", fill = TRUE) #this file need use "\t" to seperate
pca <- prcomp(gene_data[,2:51],scale = TRUE)
pcadf<- data.frame(pca$rotation)
gene_data_analysis <- cbind(samples,pcadf)
write.csv(gene_data_analysis,"gene_sample_pca.csv")
```
```{r}
#analysis the samples data into the histogram
ggplot(data = samples) +  geom_bar(mapping = aes(x = donor_vital_status), width = 0.5)
ggplot(data = samples) +  geom_bar(mapping = aes(x = disease_status_last_followup, fill = donor_vital_status), width = 0.5)
```
```{r}
#T-test
t.test(samples$donor_age_at_last_followup[samples$disease_status_last_followup=="relapse"],samples$donor_age_at_last_followup[samples$disease_status_last_followup=="complete remission"]) 
```
```{r}
#there are 50 samples in the dataset and compare them
corr_gene <- cor(gene_data[,2:51])
melted_gene<- melt(corr_gene) # the dataset will become two variable column and one value column
head(melted_gene)
p<-ggplot(melted_gene , aes(x = Var1, y = Var2)) + geom_raster(aes(fill = value)) + scale_fill_gradient2(low="navy", mid="white", high="red", midpoint=0.5) + theme( plot.title = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank())
ggplotly(p)
```
```{r}
#make a linkage map, and this can be plotted as a dendrogram
clusters <- hclust(dist(gene_data_analysis[,10:11]))
dend <- as.dendrogram(clusters)
dend <- color_branches(dend, k=3)
par(cex=0.5) # reduces font
plot(dend)
```
```{r}
#put some data into pca
write.csv(pcadf,"pcadf.csv")
plot_ly(data = gene_data_analysis, x = ~PC1, y = ~PC2, color = gene_data_analysis$icgc_donor_id )
plot_ly(data = pcadf, x = ~PC1, y = ~PC2, text = rownames(pcadf))
plot_ly(pcadf, x = ~PC1, y = ~PC2, z = ~PC3, color = ~PC4, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
 layout(scene = list(xaxis = list(title = 'PC1'),
                     yaxis = list(title = 'PC2'),
                     zaxis = list(title = 'PC3')))
```

## Python Code part

Code of creating the expression data

```{r}
import pandas as pd
import numpy as np

if __name__ == '__main__':
	data = pd.read_csv('exp_seq.BRCA-KR.tsv', sep='\t')

	n_samples = data.shape[0]
	sample_ids = data['icgc_sample_id'].unique()
	gene_ids = data['gene_id'].unique()
	col_names = np.append('gene_id', sample_ids)
	cols = col_names.shape[0]
	rows = gene_ids.shape[0]

	# np.zeros((rows, cols))
	new_data = pd.DataFrame(np.zeros((rows, cols)), columns=col_names, index=gene_ids)
	new_data = new_data.astype('str')
	# print(new_data)

	for r in range(n_samples):
	    if r % 1000 == 0:
	        print(r)
	    col = data.iloc[r]['icgc_sample_id']
	    gene_id = data.iloc[r]['gene_id']
	#     print(gene_id)
	    new_data.loc[gene_id][col] = data.iloc[r]['normalized_read_count']
	    new_data.loc[gene_id]['gene_id'] = gene_id

	new_data.to_csv('data.csv', sep='\t', index=False)
```

## Shinyapp.io Part

The sample shiny part

```{r}
library(ggplot2)
library(shiny)
library(RColorBrewer)
dataset <- read.csv('gene_sample_pca.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
headerNames=colnames(dataset)
ui <- fluidPage(
  pageWithSidebar(
    
    headerPanel("Data Explorer"),
    
    sidebarPanel(
      selectInput('x', 'X', c("None"=FALSE,headerNames),headerNames[2]),
      selectInput('y', 'Y', c("None"=FALSE,headerNames),headerNames[3]),
      selectInput('color', 'Color', c("None"=FALSE,headerNames),headerNames[5]),
      
      checkboxInput('geom_point', 'geom_point',TRUE),
      checkboxInput('geom_dotplot', 'geom_dotplot'),
      checkboxInput('geom_bar', 'geom_bar'),
      checkboxInput('geom_violin','geom_violin'),
      checkboxInput('geom_histogram','geom_histogram')
    ),
    
    mainPanel(
      plotOutput('plot')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot <- renderPlot({
    p <- ggplot(dataset, aes_string(
      x=input$x, fill=input$fill, size=input$size, color=input$color))
    #+scale_size_continuous(range=input$sizeRange)
    if (input$geom_point)
      p <- p + geom_point(aes_string(y=input$y))
  
    if (input$geom_bar)
      p <- p + geom_bar() 
    if (input$geom_dotplot)
      p <- p + geom_dotplot()
    if (input$geom_violin)
      p <- p + geom_violin(aes_string(y=input$y))
    if (input$geom_histogram)
      p <- p + geom_histogram()
    print(p)
    
  }, height=700)
  
}

# Run the application 
shinyApp(ui = ui, server = server)
```

The PCA part shiny

```{r}
library(shiny)
library(ggplot2)
library(RColorBrewer)
library(plotly)

dataset <- read.csv('gene_sample_pca.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
headerNames=colnames(dataset)

ui <- fluidPage(
  pageWithSidebar(
    
    headerPanel("PCA Explorer"),
    
    sidebarPanel(
      selectInput('x', 'X', c("None"=FALSE,headerNames),headerNames[11]),
      selectInput('y', 'Y', c("None"=FALSE,headerNames),headerNames[12]),
      selectInput('z', 'z', c("None"=FALSE,headerNames),headerNames[13])
    ),
    
    mainPanel(
      plotlyOutput('plot')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot <- renderPlotly({
    plot_ly(dataset, x = ~get(input$x), y = ~get(input$y), z = ~get(input$z), color = dataset$icgc_donor_id, colors = c('#BF382A', '#0C4B8E')) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')))
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
```

## The Website Published Part

My final project file-finalproject.Rmd has been already uploaded into Rpubs:
http://rpubs.com/vidia/FinalProjectBRCA

Also, the sample shinyapp.io published website:
https://xianghuildataexplorer.shinyapps.io/samples/

Finally, the PCA into shinyapp.io: https://xianghuildataexplorer.shinyapps.io/PCAshinyapp/
However, there is an error occured when I published. Thus the code is above.





