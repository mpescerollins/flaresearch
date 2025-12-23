import pandas as pd
from sunpy.net import Fido
from sunpy.net import attrs as a

def goes_flares_gt_M(tstart, tend):
    """
    Return a DataFrame of GOES flares with class > M (M and X) in [tstart, tend].
    Times are UTC strings as returned by HEK.
    """
    print(f"Querying GOES flares > M from {tstart} to {tend}...")
    # Query HEK for flare (FL) events with GOES class > M1.0
    res = Fido.search(
        a.Time(tstart, tend),
        a.hek.EventType("FL"),
        a.hek.FL.GOESCls > "M1.0",
        a.hek.OBS.Observatory == "GOES"
    )

    # HEK results as an Astropy table
    tab = res["hek"]
    filtered_results = tab[["event_starttime", "event_peaktime",
                               "event_endtime", "fl_goescls", "ar_noaanum"]]
    by_magnitude = sorted(filtered_results, key=lambda x: ord(x['fl_goescls'][0]) + float(x['fl_goescls'][1:]), reverse=True)
    #open a file and print the results of the query, the goes classs and start time.
    goes_str = '#Goes class, start time\n'
    for flare in by_magnitude:
        print(f"Class {flare['fl_goescls']} occurred on {flare['event_starttime']}")
        formated_time = flare['event_starttime'].strftime('%Y-%m-%dT%H:%M:%S')
        goes_str += f"{flare['fl_goescls']}, {formated_time}\n"
    #save to file
    with open("goes_flares_output.txt", "w") as f:
        print("Writing results to goes_flares_output.txt")
        f.write(goes_str)
    

if __name__ == "__main__":
    # Example usage
    #
    # Get all GOES flares > M between March 1, 2024 and May 1, 2024
   

    #res = Fido.search(a.Time("2011/08/09 07:23:56", "2011/08/09 12:40:29"),
    #                  a.hek.EventType("FL"))
    #print(res)
    # 
    #df = goes_flares_gt_M("2024-03-01", "2024-05-01")

    goes_flares_gt_M("2024-03-01", "2024-04-03")
    # Your "list of flare times" (begin times) for flares > M:
    #flare_start_times = df["begin_time"].tolist() if "begin_time" in df.columns else []
    #print(flare_start_times)