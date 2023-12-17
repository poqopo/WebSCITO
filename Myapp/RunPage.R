header <-div(class="sticky top-0 flex place-content-between w-full h-[100x] border-b-2 px-10 border-black bg-white z-50",
             img(class='my-auto w-1/4 max-w-[250px] max-h-[80px] object-cover text-[30px] font-extrabold text-start', src='logo.png'),
             div(class='my-auto w-1/2 flex place-content-between text-[20px] text-center font-extrabold',
                 a(href=route_link("/"),"Introduction"),
                 a(href=route_link("tutorial"),"Tutorial"),
                 a(href=route_link("help"),"Help"),
                 a(href=route_link("contact"), "Contact"),
                 a(href=route_link("run"), "RUN")))


hashtag_loading <-  div(class = "mx-auto my-auto flex place-content-center",
                        p(class = "font-bold text-[20px] px-10 my-auto w-fit", "Making Hashtag Object..."),
                        img(class='w-[50px] h-[50px]', src='loading.gif'))
multimodal_loading <-  div(class = "mx-auto my-auto flex place-content-center",
                           p(class = "font-bold text-[20px] px-10 my-auto w-fit", "Making multimodal Seurat Object..."),
                           img(class='w-[50px] h-[50px]', src='loading.gif'))
resolve_loading <-  div(class = "mx-auto my-auto flex place-content-center",
                        p(class = "font-bold text-[20px] px-10 my-auto", "Making Resolved Object..."),
                        img(class='w-[50px] h-[50px]', src='loading.gif'))
integration_loading <-  div(class = "mx-auto my-auto flex place-content-center",
                            p(class = "font-bold text-[20px] px-10 my-auto", "Batch Specific Normalize and Integrating Input..."),
                            img(class='w-[50px] h-[50px]', src='loading.gif'))


run_page = div(
  header,
  div(class="w-2/3 max-w-[1200px] p-5 mx-auto",
      h3(class = "py-10 font-extrabold text-[30px] text-center", "Start your own research"),
      div(class ='w-full flex place-content-between',
          div(class ='m-auto w-4/5 flex place-content-evenly gap-x-2',
              fileInput(inputId="raw", accept = c(".h5", ".rds") ,label='filtered_h5_bc_matrix file input', width = "25%"),
              numericInput('abNum', label = "Antibody number", value = 1, width = "20%"),
              selectInput("batchNum", "Batch Number", choices = c(1:10), width = '20%'),
              div(class= "text-[20px]",checkboxInput("multimodal", label = "Multimodal data", value = FALSE)),
          ),
          actionButton(class="mx-10 w-full max-w-[200px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-full",
                       inputId = "start", label = "Start")
      ),
      div(id = "hashtagloading", hashtag_loading),
      uiOutput("basicInfo"),
      div(id = "multimodalloading", multimodal_loading),
      uiOutput("ridgeplotUI"),
      div(id = "resolveloading", resolve_loading),
      div(id = "integrationloading", integration_loading),
      uiOutput("pcaumap"),
      uiOutput("dotplotUI"),
      uiOutput("featureplotUI")
  )
)
