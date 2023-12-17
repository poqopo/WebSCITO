header <-div(class="sticky top-0 flex place-content-between w-full h-[100x] border-b-2 px-10 border-black bg-white z-50",
             img(class='my-auto w-1/4 max-w-[250px] max-h-[80px] object-cover text-[30px] font-extrabold text-start', src='logo.png'),
             div(class='my-auto w-1/2 flex place-content-between text-[20px] text-center font-extrabold',
                 a(href=route_link("/"),"Introduction"),
                 a(href=route_link("tutorial"),"Tutorial"),
                 a(href=route_link("help"),"Help"),
                 a(href=route_link("contact"), "Contact"),
                 a(href=route_link("run"), "RUN")))


contact_page = div(
  header,
  div(class="w-2/3 max-w-[1200px] m-auto ",
      h1(class="text-[30px] font-extrabold py-5" ,"Contact US"),
      h3(class="text-[24px] font-bold leading-normal py-5", "If you have any questions, Please shoot us an email:)"),
      h3(class="text-[18px] font-semibold leading-normal pt-5", "Main Contributor :"),
      h3(class="text-[18px] leading-normal",
         span("Sangjun Park(Undergraduate, YONSEI University) :"),
         span(class = "underline underline-offset-1", "bio.poqopo.gmail.com")),
      h3(class="text-[18px] font-semibold leading-normal pt-5", "Corresponding author :"),
      h3(class="text-[18px] leading-normal",
         span("Byungjin Hwang(Assistant Professor, Yonsei University, College of Medicine) :"),
         span(class = "underline underline-offset-1", "	bjhwang113@yuhs.ac")),
      p(class = "text-[18px] pt-5 font-semibold", "Visit our website for more information!"),
      a(class = "py-5 text-[18px] underline", href="https://sites.google.com/view/bhwanglabyonsei/home",
        "https://sites.google.com/view/bhwanglabyonsei/home")
    
  )
)
