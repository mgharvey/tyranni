pkgconfig::set_config("drake::strings_in_dots" = "literals")

my_plan <- drake_plan(
  phy = GetTree(),
  sisters = sis_get_sisters(phy, ncores=1),
  traits = GetData(),
  cleaned = CleanData(phy, traits),
  phy_cleaned = cleaned$phy,
  traits_cleaned = cleaned$traits,
  sister_comparisons_25 = DoAllTraits(sisters, traits_cleaned, phy_cleaned, cutoff=0.25),
  sister_comparisons_50 = DoAllTraits(sisters, traits_cleaned, phy_cleaned, cutoff=0.5),
  sister_comparisons_75 = DoAllTraits(sisters, traits_cleaned, phy_cleaned, cutoff=0.75),
  sister_combined = CombineResults(sister_comparisons_25, sister_comparisons_50, sister_comparisons_75),
  saved = write.csv(sister_combined, file_out("Sister.csv")),
  heatmap = DoHeatmap(sister_combined),
  iterated = sisters::sis_iterate(traits_cleaned[,"species.richness"], nsteps=101, phy=phy_cleaned, sisters=sisters),
  plotpairsresult = PlotPairs(iterated),
  save_iterated = write.csv(iterated, file_out("Sensitivity.csv"))
)
