# sediment_screening_tool
this is a draft as of 8/25/2023
 NOT THE FINAL VERSION - NEEDS FURTHER TESTING TO BE SURE EVERYTHING IS CORRECT

Basic sediment metrics
1.	Marble Canyon sand mass balance – computed from a modified form of the Wright et al. (2010) sand routing model that can run on monthly volumes. Uses average annual Paria sediment input (using most-likely value, not lower bound). Does not include HFEs (because HFEs are not in annual/monthly data). Mass balance is computed on annual (calendar year) timestep.
2. LTEMP HFE metric. HFE implementation probability per year for sediment year (july 1st to june 30th). Probability that mass balance is sufficient to support a 60-hr HFE and the following constraints are met:
  - Assume Powell reservoir elevation must be at least 3550' for fall implementation and 3525' for spring, because HFE's have not been implemented below these elevations
  - All HFEs are assumed to be 40,000 cfs magnitude and at least 60 hr duration
  - The annual HFE implementation probability will be used in Brad Butterfield's vegetation modeling
3. SEIS metric - probably not included in the screening tool but computed in the code anyways. Similar to above but using the revised 1-yr sediment accounting window rather than LTEMP windows.
  - fall implementation same as above
  - spring implementation months (apr, may, jun) computed with 1-yr sediment accounting window and 3525' reservoir constraint

Procedure:
1.	Read output from CRSS
-	Lake Powell end of month elevations.
- Glen Canyon Dam monthly release volume.
3.	Run monthly version of sand routing model to compute metric 1 (Marble Canyon sand mass balance)
   - Assume a fixed bed composition equal to the median bed composition computed by the Wright et al. (2010) model for 2002-2023. (Without including HFE’s explicitly, the bed condition cannot be known because HFE’s coarsen the bed, therefore it is better to      
     assume a fixed “typical” bed condition.)
   - Calculate daily high and low flows based on LTEMP constraints
   - Use relation from Wright et al. (2010) to compute the monthly sand export from Marble Canyon based on the bed condition specified above, and assuming 50% of the time at high flow and 50% at low flow.
   - Use the average monthly Paria sediment load to compute the cumulative monthly sand mass balance for Marble Canyon (assuming additional tributary sand inputs equal to 10% of Paria). The average monthly Paria sediment loads will be hard-coded into the script       to avoid reading an external file.
   - Output for sand mass balance metric: Annual mass balance (calendar year)
   - Possible ways to turn this into a performance metric:
      - average annual mass balance across all years and all traces
      - percent of years below threshold mass balance (294 Tmt, which is equal to a 60-hr HFE) across all years and all traces
  5.	Compute LTEMP HFE metric
   - For each year in the 30-year trace, find the sand export through the end of November (computed in step 2) plus the sand export resulting from a 60-hr 40,000 cfs HFE (294 Tmt based on Wright et al. model using bed condition from step 2)
   - Based on historical Paria inputs, find the percentage of years where sediment supply from July 1 to November 1 exceeds the sand export (including HFE). This represents the probability of a sediment-triggered HFE for that year. The 26 values for Paria             sediment supply are read from external file. Similar for the spring implementation window.
   - If the October end-of-month Lake Powell elevation is below 3550 ft, assume that an HFE cannot be implemented, and therefore set the HFE implementation probability to 0. Otherwise, the HFE implementation probability is the same as the probability of a               sediment-triggered HFE for that fall. Similar for spring, if march end-of-month powell elevation is below 3525 ft assume HFE cannot be implemented.
   - Annual probability is the probability that a fall HFE AND/OR a spring HFE is implemented
   - Potential performance metric: average annual probablity of HFE implementation across all years and all traces
  6.  Compute SEIS HFE metric (optional) - similar to above but 1-yr window and four possible implementation months, annual probability is max probability for any of the four implementation months (probabilities not independent so no AND/OR)
