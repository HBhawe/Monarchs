table = read.csv("migrant_vs_resident_Final_stats.csv")

test1 = table[order(-table$West_vs_resident_Fst), ]
westvresident = test1[1:236, ]
View(westvresident)

test2 = table[order(-table$East_vs_resident_Fst), ]
eastvresident = test2[1:236, ]

test3 = table[order(-table$migrant_vs_resident_Fst), ]
migrantvresident = test3[1:236, ]

top2 = test3[1:472, ]
top2_migrantvresident = write.csv(top2, "top2migrantvresident.csv")

top5 = test3[1:1180,]
top5_migrantvresident = write.csv(top5, "top5migrantvresident.csv")

florida_west = table[order(-table$West_vs_Florida_FST), ]
florida_east = table[order(-table$East_vs_Florida_FST), ]


write.csv(westvresident, "West_vs_Resident.csv")
write.csv(eastvresident, "East_vs_Resident.csv")
write.csv(migrantvresident, "migrant_vs_Resident.csv")
write.csv(florida_west, "westvflorida.csv")
write.csv(florida_east, "eastvflorida.csv")
