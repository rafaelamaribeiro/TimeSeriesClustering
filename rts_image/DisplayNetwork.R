# R Netplan

library(shiny)
library(leaflet)
library(rgdal)
library(leaflet.extras)

CreateSLDF <- function(netData) {
  sldf <- lapply(X = 1:nrow(netData), FUN = function(x) {
    l <- Line(coords = rbind(as.numeric(netData[x,c('Lon.x','Lat.x')]),
                             as.numeric(netData[x,c('Lon.y','Lat.y')]))
    )
    
    ls <- Lines(slinelist = l, ID = x)
  })
  
  
  sldf <- SpatialLines(LinesList = sldf)
  
  nl <- nrow(netData)
  sldf <- SpatialLinesDataFrame(sl = sldf, data = netData, match.ID = F)
  
  return(sldf)
}

ui <- fluidPage(
  leafletOutput("NetViewR",width = "100%",height = "800px")
)

server <- function(input, output, session) {

  output$NetViewR <- renderLeaflet({
  
    NetMap <- leaflet(data = bus)
    #NetMap <- NetMap %>% addProviderTiles(providers$Esri.WorldStreetMap,
    #                                      options = providerTileOptions(noWrap = TRUE))
    NetMap <- NetMap %>% addCircles(~Lon, ~Lat, color="black", popup=~paste0('<b> ID: ',ID,'</b>'), weight = 5, radius=20, stroke = TRUE, fillOpacity = 1.0)
    NetMap <- NetMap %>% addPolylines(data = cirSLDF, color="black", weight = 1, popup=~paste0('<b> ID FROM: ',ID.FROM,'</b><\br><b> ID TO:',ID.TO,'</b>'))
    # NetMap
    NetMap <- NetMap %>% setMapWidgetStyle(list(background= "white"))}

  )
}

bus <- read.csv("C:\\Users\\frog-\\OneDrive\\Rafaela\\TimeSeriesClustering\\Nodes_RTS-GMLC.csv")
cir <- read.csv("C:\\Users\\frog-\\OneDrive\\Rafaela\\TimeSeriesClustering\\Lines_RTS-GMLC.csv")

cirXY <- merge(x = cir  , y = bus, by.y="ID", by.x="ID.FROM")
cirXY <- merge(x = cirXY, y = bus, by.y="ID", by.x="ID.TO"  )

cirSLDF <- CreateSLDF(cirXY)

plot(cirSLDF)

shinyApp(ui, server)
