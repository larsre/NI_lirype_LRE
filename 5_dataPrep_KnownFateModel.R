

###########################################################################################
#### DATA FROM RADIO TELEMETRY PROJECT; 


library(RSQLite)
library(lubridate)
library(tidyverse)

con <- dbConnect(drv=RSQLite::SQLite(), dbname="P:/12179000_Lirypetelemetri_i_Lierne/DATA Lierne/Database/Lierne_WP_Project_Data_ReSt.db")
Occur <- as_tibble(dbReadTable(con, "Marked_occurences"))
Cap <- as_tibble(dbReadTable(con, "Marked_CaptureDetails"))

dbDisconnect(con)


####################################
### Preparing for analysis; 

d_m <- Cap %>% select(organismID, CaptureAge, Sex_fieldbased, Sex_dna, CaptureArea) %>%
        filter(CaptureArea!="Fiskløysdalen") %>%
        left_join(., Occur) %>%
        select(organismID, CaptureAge, Sex_dna, Sex_fieldbased, occurenceID, Date, Status, Comment, RingNR, CaptureArea) %>%
        mutate(sex=if_else(!is.na(Sex_dna), Sex_dna, Sex_fieldbased)) 
write.csv(d_m, "demographic_data/ObsData.csv", row.names = FALSE)
rm(d_m)


############################################################################




d_m <- as_tibble(read.csv("demographic_data/ObsData.csv", header=T, sep=",")) %>%
        mutate(Date=lubridate::dmy(Date)) %>%
        mutate(YearDay=yday(Date), Year=year(Date), JD=julian(Date, origin=as.Date("2015-08-01"))) %>%
        mutate(yday2=if_else(YearDay>yday(ymd("2015-08-01")), YearDay-yday(ymd("2015-08-01")), YearDay+(365-yday(ymd("2015-08-01"))))) %>%
        mutate(YearPeriod=if_else(YearDay<yday(ymd("2015-08-01")), Year-1, Year)) %>%
        mutate(TimePeriod=cut(yday2, breaks=c(0, 180, 365), labels=c(1,2)), TimePeriod=as.numeric(TimePeriod))


### Deaths
temp_dead <- d_m %>% filter(between(Status, 3,6)) %>%
            group_by(organismID) %>%
            filter(Date==min(Date)) %>%
            slice(1) %>%
            ungroup() %>%
            select(RingNR, sex, Date, YearPeriod, Year, TimePeriod, yday2, Status) 

### Captures; 

temp_cap <- d_m %>% filter(Status==1) %>%
           select(RingNR, sex, Date, YearPeriod, Year, TimePeriod, yday2)
         

### Live

temp_live <- d_m %>% filter(between(Status, 1,2)) %>%
         select(RingNR, sex, Date, YearPeriod, Year, TimePeriod, yday2, Status) %>%
         bind_rows(., temp_dead) %>%
         group_by(RingNR, YearPeriod, TimePeriod) %>%
         summarize(Antall=n()) %>%
         ungroup() %>%
         mutate(Y_P=YearPeriod+(0.1*TimePeriod)) %>%
         group_by(RingNR) %>%
         mutate(LastPeriod=max(Y_P), Excl=if_else(Y_P==LastPeriod, 0,1)) %>%
         ungroup() %>%
         filter(Excl==1) %>%
         mutate(Incl=1, Faith=1) %>%
         select(RingNR, YearPeriod, TimePeriod, Incl, Faith)




CaptureHistory <- temp_dead %>% select(RingNR, YearPeriod, TimePeriod) %>%
                 mutate(Incl=1, Faith=0) %>%
                 bind_rows(., temp_live) %>%
                 mutate(Y_P=YearPeriod+(0.1*TimePeriod)) %>%
                filter(between(Y_P, 2015.1, 2019.9))

CH_m <- CaptureHistory %>% group_by(YearPeriod, TimePeriod) %>%
        summarize(AtRisk=sum(Incl), Survivors=sum(Faith)) %>%
        ungroup()

filter(CH_m, TimePeriod==1)
filter(CH_m, TimePeriod==2)




####################################################################

test <- d_m %>% 
        filter(Status!=7) %>%
        group_by(RingNR) %>%
         mutate(LastPeriod=max(Date),  Excl=if_else(Date==LastPeriod, 0,1)) %>%
         ungroup() %>%
         filter(Excl==0, Status==2) %>%
        select(RingNR, Date, Status)

########################################################################

tempo <- d_m %>% filter(between(Status, 1,2)) %>%
  select(RingNR, sex, Date, YearPeriod, Year, TimePeriod, yday2, Status) %>%
  bind_rows(., temp_dead) %>% 
  group_by(RingNR, YearPeriod, TimePeriod) %>%
  summarize(Antall=sum(Status), Antall2=n()) %>%
  complete(TimePeriod=1:3, nesting(YearPeriod, RingNR), fill=list(Antall=0, Antall2=0)) %>%
  ungroup() %>%
  mutate(ID=row_number()) 

test <- tempo %>% group_by(RingNR) %>%
        mutate(Først=first(ID, Antall>0)) 



            
