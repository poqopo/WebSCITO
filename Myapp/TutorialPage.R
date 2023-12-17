
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
                           p(class = "font-bold text-[20px] px-10 my-auto w-fit", "Making Multimodal Object..."),
                           img(class='w-[50px] h-[50px]', src='loading.gif'))
resolve_loading <-  div(class = "mx-auto my-auto flex place-content-center",
                        p(class = "font-bold text-[20px] px-10 my-auto", "Making Resolved Object..."),
                        img(class='w-[50px] h-[50px]', src='loading.gif'))
integration_loading <-  div(class = "mx-auto my-auto flex place-content-center",
                            p(class = "font-bold text-[20px] px-10 my-auto", "Split and Integrating Input..."),
                            img(class='w-[50px] h-[50px]', src='loading.gif'))



tutorial_page <- div(
  header,
  div(class="w-2/3 max-w-[1200px] p-5 mx-auto",
      div(
        h3(class = "py-10 font-extrabold text-[30px] text-start", "Select Test Samples"),
        div(class ='my-10 w-full flex place-content-evenly',
            selectInput("Samples", "Samples:",
                        list(`Unimodal Data` = list("100K_PBMC_10batch" = "/data/project/scitoseq/script/Final/Sample/100K_28AB_10Pool/",
                                                    "200K_PBMC_10batch" = "/data/project/scitoseq/script/Final/Sample/200K_28AB_10Pool/"),
                             `Multimodal Data` = list("PBMC_4Batch" = "/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/"))),
            actionButton(class="mx-10 w-full max-w-[200px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-full",
                         inputId = "tutorial_start", label = "Start")
        ),
        div(id = "shashtagloading", hashtag_loading),
        div(id = "smultimodalloading", multimodal_loading),
        uiOutput("sbasicInfo"),
        uiOutput("sridgeplotUI"),
        div(id = "sresolveloading", resolve_loading),
        div(id = "sintegrationloading", integration_loading),
        uiOutput("spcaumap"),
        uiOutput("sdotplotUI"),
        uiOutput("sfeatureplotUI")
      ),
  )
)