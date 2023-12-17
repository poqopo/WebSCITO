header <-div(class="sticky top-0 flex place-content-between w-full h-[100x] border-b-2 px-10 border-black bg-white z-50",
             img(class='my-auto w-1/4 max-w-[250px] max-h-[80px] object-cover text-[30px] font-extrabold text-start', src='logo.png'),
             div(class='my-auto w-1/2 flex place-content-between text-[20px] text-center font-extrabold',
                 a(href=route_link("/"),"Introduction"),
                 a(href=route_link("tutorial"),"Tutorial"),
                 a(href=route_link("help"),"Help"),
                 a(href=route_link("contact"), "Contact"),
                 a(href=route_link("run"), "RUN")))

help_page = div(
  header,
  div(class="w-2/3 max-w-[1200px] mx-auto",
      div(id = "Explain",
          h1(class = "font-extrabold text-[36px] my-10",
             "WORKFLOW"),
          img(src='WORKFLOW.png', class = "h-[400px] mx-auto object-cover", alt = "Loading..."),
          h2(class = "font-bold text-[30px] my-5", "Input File"),
          
          
          p(class= "text-[18px] mt-10",
            "The main objective of WebSCITO is to analyze scitoseq data. This method utilizes a 
            combination of split-pool techniques and droplet single-cell combinatorial sequencing 
            (dsc-seq), merging antibody barcodes with pool barcodes to enable the analysis of 100K-1M 
            cells per reaction. This approach marks a notable advancement in single-cell analysis, 
            aiming to deliver results with enhanced precision and efficiency."),
          p(class= "text-[18px] mt-5 mb-10",
            "WebSCITO also assists users in analyzing scito-seq data by addressing issues like batch 
            effects and integrating multimodal data. For guidance on preparing samples, please visit 4
            the provided link."),
          a(href = "https://www.nature.com/articles/s41592-021-01222-3", class = "text-[16px] font-bold my-1",
            "SCITO-seq: single-cell combinatorial indexed cytometry sequencing"),
          
          
          h2(class = "font-bold text-[24px] my-10", "Filtered Cell Matrix(.h5 file)"),
          p(class= "text-[18px]",
            "On the Run page, you can begin the analysis using your own data. Start with the initial 
            input, which is the filtered_bc_matrix, an output from cellRanger. Once you see the 'Upload 
            Complete' message, proceed to provide additional information such as antibody and batch 
            numbers. To initiate WebSCITO, click the 'Start Run' button. If your data is multimodal, 
            make sure to check the corresponding checkbox."),
          h2(class = "font-bold text-[24px] my-10", "Antibody Number and Batch Number"),
          p(class= "text-[18px]",
            "The Antibody Number corresponds to the total count of antibodies used, representing the 
            overall number of unique antibodies in your analysis. The Batch Number indicates how the 
            samples are divided into batches for processing. Specifically, the Antibody Number is 
            crucial for combining pool-specific antibodies, while the Batch Number is essential for 
            normalizing data within each batch."),
          p(class= "text-[18px] mt-5",
            "Ensure that you input the correct Antibody Number and Batch Number. For example, if you are 
            using 28 antibodies distributed across 10 batches, the total number of rows in your h5 file 
            should be 280 (28 * 10). In this case, the Antibody Number would be 28, and the Batch Number 
            would be 10. Make sure to accurately specify these values."),
          img(src='helpinput.png', class = "w-full mx-auto object-cover my-10 ", alt = "Loading..."),
          h2(class = "font-bold text-[24px] my-10", "Batch Metadata (For Multiomics)"),
          p(class= "text-[18px]",
            "In WebSCITO v1, during multimodal data analysis, RNA barcodes are utilized to categorize 
            cells into distinct batches. For this reason, users are required to provide batch information
            related to the cells. In contrast, for unimodal data analysis, WebSCITO will automatically 
            assign each batch based on antibody expression. Users do not need to manually specify batch 
            information in unimodal analysis"),
          h2(class = "font-bold text-[24px] my-10", "Metadata File (Optional)"),
          p(class= "text-[18px] mb-5",
            "Users can include additional metadata for their data at this stage. This can involve adding 
            cell metadata, such as cell type annotations. It's important to ensure that the row number in 
            the metadata file matches the number of input cells. For reference, an example of a metadata 
            file is provided below."),
          dataTableOutput("helpmeta"),
          p(class= "text-[18px] mt-5",
            "Default format of metadata file is 'csv'"),
          h2(class = "font-bold text-[24px] my-10", "Antibody Name File (Optional)"),
          p(class= "text-[18px]",
            "The Antibody Name File is an optional file. During the integration of each batch, WebSCITO 
            uses pool-specific antibodies, resulting in resolved antibody names that may slightly differ 
            from the actual names (e.g., P1_CD45). Providing the actual names ensures a more seamless 
            representation in the visualization. It's crucial to note that the number of antibody names 
            in this file should match the antibody number specified in the first step of the process."),
          h1(class = "font-extrabold text-[36px] my-10",
             "Step"),
          h2(class = "font-bold text-[30px] my-10", "1. Make Hashtag Object"),
          p(class= "text-[18px]",
            "In the initial stage of the WebSCITO process, the hashtag object is created. During this 
            step, cells are allocated to specific batches based on their antibody expressions. It's 
            important to note that in the case of multiplets, a single droplet can be assigned to 
            multiple batches. To provide users with insights into the batch resolution results and 
            efficiency, visual representations such as pie charts and ridge plots are available."),
          p(class= "text-[18px] mt-5",
            "However The process differs in multimodal analysis. In multimodal analysis, cells are 
            allocated based on RNA expression. Therefore, the creation of the hashtag object step is 
            skipped in multimodal analysis."),
          h2(class = "font-bold text-[30px] my-10", "2. Make Resolved Object"),
          p(class= "text-[18px]",
            "In this step, WebSCITO generates a resolved object utilizing the hashtag object's batch 
            allocation. The resolved object is created by resolving inputs based on their respective 
            batches. It's important to note that due to the assignment of multiplets to multiple batches,
            the resolved object may have more cells than the original hashtag object. Additionally, 
            during this process, split antibodies are combined together."),
          h2(class = "font-bold text-[30px] mt-10", "3. Quality Control and Normalization"),
          p(class= "text-[18px] pt-10",
            "After the creation of the resolved object, WebSCITO advances to the normalization step. 
            To address batch effects, WebSCITO divides the data by batch and independently conducts 
            normalization for each batch. Subsequently, the normalized results are integrated. WebSCITO 
            employs rpca normalization methods to integrate results."),
          p(class= "text-[18px] pt-5 pb-10",
            "In the context of multimodal analysis, WebSCITO utilizes the WNN approach to integrate both 
            RNA and ADT data. For more detailed information regarding the process of data splitting and 
            integration, please consult the provided link."),
          a(href = "https://satijalab.org/seurat/articles/seurat5_integration", class = "text-[16px] font-bold my-10",
            "Integrative analysis in Seurat v5"),
          h2(class = "font-bold text-[30px] my-10", "4. Visualization"),
          p(class= "text-[20px] font-bold py-10", "1) Pie Chart and Ridge Plot"),
          p(class= "text-[18px]",
            "This chart provides an overview of hashtag object resolution through batch processing. The 
            first two pie charts illustrate the distribution of batch types and the corresponding cell 
            counts. These charts offer insights into the composition of the resolution process. 
            Additionally, the ridgeplot visually represents the efficiency of the resolution across 
            batches."),
          div(class = "w-full flex place-content-evenly",
              img(src='batch percentage.png', class = "h-[300px] object-cover my-10 ", alt = "Loading..."),
              img(src='singlet percentage.png', class = "h-[300px] object-cover my-10 ", alt = "Loading...")
          ),
          img(src='ridgeplot.png', class = "w-full object-cover my-10 ", alt = "Loading..."),
          p(class= "text-[20px] font-bold py-10", "2) PCA and UMAP"),
          p(class= "text-[18px]",
            "Following the normalization process, WebSCITO displays a PCA elbow plot and UMAP clustering 
            visualization of the results. The PCA elbow plot provides a clear depiction of the 
            effectiveness of Principal Components (PCs) definition. Simultaneously, the UMAP clustering 
            visualization showcases the grouping of cells, offering insights into the inherent patterns 
            within the data. Notably, users have the capability to switch between different groups in 
            the UMAP visualization, a feature provided in additional metadata. "),
          img(src='pca_umap.png', class = "w-full object-cover my-10 ", alt = "Loading..."),
          p(class= "text-[20px] font-bold py-10", "3) Dotplot"),
          p(class= "text-[18px]",
            "WebSCITO shows a dot plot that illustrates individual antibody expressions for each 
            group. Similar to the UMAP visualization, users have the capability to switch between 
            different groups (e.g., Batch, Samples) to examine expression patterns. Moreover, users 
            can selectively choose antibodies of interest to explore specific features and groups. 
            In the context of multimodal data, WebSCITO also provides an RNA dot plot, offering a 
            comprehensive view of RNA expression patterns with antibodies."),
          img(src='Dotplot.png', class = "w-full object-cover my-10 ", alt = "Loading..."),
          p(class= "text-[20px] font-bold py-10", "4) Featureplot"),
          p(class= "text-[18px]",
            "WebSCITO provides users with the capability to investigate individual antibody expressions 
            within the UMAP using feature plots. Similar to other visualizations, users can choose 
            specific antibodies for observation. In the case of multimodal data, WebSCITO enhances this 
            functionality by incorporating feature plots for RNA data, enabling users to visualize 
            RNA expressions. This comprehensive approach allows users to compare and contrast the 
            two plots, facilitating a robust understanding of both antibody and RNA expression 
            patterns within the context of the UMAP."),
          img(src='FeaturePlot.png', class = "w-full object-cover my-10 ", alt = "Loading...")
      )
  )
)