from python_cc_reader.inclusion_removal.count_total_inclusions import total_inclusion_count


if __name__ == "__main__":
    internal_count, all_count = total_inclusion_count()
    print("internal nonfwd includes:", internal_count, "all nonfwd includes:", all_count)
