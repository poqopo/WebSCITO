header <-div(class="sticky top-0 flex place-content-between w-full h-[100x] border-b-2 px-10 border-black bg-white z-50",
             img(class='my-auto w-1/4 max-w-[250px] max-h-[80px] object-cover text-[30px] font-extrabold text-start', src='logo.png'),
             div(class='my-auto w-1/2 flex place-content-between text-[20px] text-center font-extrabold',
                 a(href=route_link("/"),"Introduction"),
                 a(href=route_link("tutorial"),"Tutorial"),
                 a(href=route_link("help"),"Help"),
                 a(href=route_link("contact"), "Contact"),
                 a(href=route_link("run"), "RUN")))

root_page = div(
  header,
  div(class="relative w-full h-[300px]",
      img(src='background3.jpg', class = "w-full h-[300px] object-cover", alt = "Loading..."),
      div(class="absolute inset-0 flex flex-col items-center justify-center drop-shadow-2xl",
          p(class="text-white text-[50px] font-extrabold", "WebSCITO"),
          p(class="text-white text-[30px] font-extrabold", "Ultra high-throughput analysis of single cell multiomics 
            with combinatorial indexing."),
  ),
  div(class="w-2/3 max-w-[1200px] m-auto ",
      h2(class="text-[30px] font-bold py-10", "Introduction to WebSCITO: Analyzer for SCITO-seq"),
      h3(class="text-[20px] leading-normal","In the rapidly evolving field of single-cell analysis, various 
         methods strive to enhance accuracy and efficiency. One standout approach is the SCITO-seq, a pioneering
         method that combines split-pool methods and droplet single-cell combinatorial sequencing (dsc-seq). By 
         combining antibody barcodes with pool barcodes, SCITO-seq enables combinatorially indexed dsc-seq of 
         antibodies, allowing the analysis of 100K-1M cells per reaction—a significant improvement over 
         conventional single-cell methods. This groundbreaking method not only traces antibodies to specific 
         pools but also effectively resolves doublets into singlets."),
      h3(class="text-[20px] leading-normal py-10","In response to the growing demand for robust SCITO-seq data
         interpretation tools, we introduce", span(class = "font-bold","WebSCITO: Ultra high-throughput analysis 
          of single cell multiomics with combinatorial indexing."), "WebSCITO seamlessly integrates with the 
         strengths of SCITO-seq."),
      h3(class="text-[20px] leading-normal","A standout feature of WebSCITO lies in its capacity to effectively 
         separate cells from different pools within doublets, utilizing ab+pool barcodes while efficiently 
         handling 100K-1M cells. This capability significantly enhances the precision and efficiency of cell 
         classification, marking a critical advancement in the field. Additionally, WebSCITO addresses batch 
         effects associated with each pool, ensuring a more thorough comprehension of the underlying biological 
         information."),
      h3(class="text-[20px] leading-normal my-10","Moreover, WebSCITO renders the seamless integration of ADT 
      and RNA data, enabling researchers to analyze multiome data and explore chosen features with a simple 
      one-click operation. This unique capability empowers scientists to gain a comprehensive understanding of 
      complex cellular diversity."),
      h3(class="text-[20px] leading-normal","In essence, WebSCITO not only conducts high-throughput analysis but
         also provides diverse visual representations to demonstrate its efficacy and facilitate subsequent 
         downstream analysis in single-cell protein and multimodal profiling. By providing large-scale 
         comprehensive single-cell multiomics analysis, WebSCITO stands at the forefront of cutting-edge 
         research in single-cell analysis."),
      h2(class="text-[30px] font-bold py-10", "Key Features of : WebSCITO"),
      div(class = "w-full flex place-content-between gap-x-10 py-10",
          div(
            p(class="text-[24px] font-bold text-center", "Multiplet resolution"),
            img(src='Doublet.png', class = "w-full max-w-[300px] my-auto", alt = "Loading..."),
          ),
          div(
            p(class="text-[24px] font-bold text-center", "Batch Effect Correction"),
            img(src='Batch.png', class = "w-full max-w-[300px] my-auto", alt = "Loading..."),            
          ),
          div(
            p(class="text-[24px] font-bold text-center", "Multimodal Integration"),
            img(src='multi.png', class = "w-full max-w-[600px] my-auto", alt = "Loading..."),
          ))
  ))
)