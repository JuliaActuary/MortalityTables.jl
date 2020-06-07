using MortalityTables, Plots


tables = MortalityTables.tables()
cso_2001 = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
cso_2017 = tables["2017 Loaded CSO Preferred Structure Nonsmoker Super Preferred Male ANB"]

issue_ages = 18:80
durations = 1:40

a = [cso_2001.select[issue_age][issue_age + duration - 1] for issue_age in issue_ages,duration in durations]
b = [cso_2017.select[issue_age][issue_age + duration - 1] for issue_age in issue_ages,duration in durations]


# compute the relative rates with the element-wise division ("brodcasting" in Julia)
function rel_diff(a, b, issue_age,duration)
        att_age = issue_age + duration - 1
        return a[issue_age][att_age] / b[issue_age][att_age]
end


diff = [rel_diff(cso_2017.select,cso_2001.select,ia,dur) for ia in issue_ages, dur in durations]
x = contour(durations,
        issue_ages,
        diff,
        xlabel="duration",ylabel="issue ages",
        title="Relative difference between 2017 and 2001 CSO \n M PFN",
        fill=true
        )

savefig(x,"contour.png")
