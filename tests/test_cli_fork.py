"""Tests for cli._seed_fork_source() (the --fork-from mechanism).

Only the pure filesystem-copying helper is unit-tested here. main()'s own
CLI orchestration (argument parsing, --fork-from validation messages,
wiring into single-/two-list mode) is not unit-tested elsewhere in this
project either — consistent with the existing convention of testing
library-level functions directly and verifying CLI-entry-point behaviour
manually/end-to-end.
"""

import os
import json

import pytest

from enhydra.cli import _seed_fork_source, _FORK_SEEDABLE_STAGE_DIRS
from enhydra.exceptions import EnhydraIOError
from enhydra.stage_markers import write_stage_marker, MARKER_FILENAME


def _touch(path):
    with open(path, "w") as fh:
        fh.write("x")


class TestSeedForkSource:

    def test_copies_all_present_stage_directories(self, tmp_path):
        src = tmp_path / "old_run"
        dst = tmp_path / "new_run"
        src.mkdir()

        for stage in ("length_filter", "group_filter", "alignment"):
            d = src / stage
            d.mkdir()
            _touch(d / "OG0001")

        _seed_fork_source(str(src), str(dst))

        for stage in ("length_filter", "group_filter", "alignment"):
            assert (dst / stage / "OG0001").is_file()

    def test_missing_stage_directories_silently_skipped(self, tmp_path):
        """A source that only completed a few early stages (e.g. an
        interrupted run) must seed only what actually exists, without
        raising for the stages that never ran."""
        src = tmp_path / "old_run"
        dst = tmp_path / "new_run"
        src.mkdir()
        (src / "length_filter").mkdir()
        _touch(src / "length_filter" / "OG0001")
        # group_filter, alignment, etc. never ran in this source.

        _seed_fork_source(str(src), str(dst))

        assert (dst / "length_filter" / "OG0001").is_file()
        assert not (dst / "group_filter").exists()
        assert not (dst / "alignment").exists()

    def test_stage_markers_survive_the_copy(self, tmp_path):
        """The .stage_marker.json inside each copied stage directory must
        be preserved intact, since _run_single_list()'s own marker check
        depends on reading it from the destination after seeding."""
        src = tmp_path / "old_run"
        dst = tmp_path / "new_run"
        src.mkdir()
        length_filter_dir = src / "length_filter"
        write_stage_marker(str(length_filter_dir), "length_filter",
                          {"length_filter_sd": 2.0})
        _touch(length_filter_dir / "OG0001")

        _seed_fork_source(str(src), str(dst))

        marker_path = dst / "length_filter" / MARKER_FILENAME
        assert marker_path.is_file()
        with open(marker_path) as fh:
            marker = json.load(fh)
        assert marker["parameters"] == {"length_filter_sd": 2.0}

    def test_only_known_stage_directories_are_copied(self, tmp_path):
        """A stray subdirectory not in _FORK_SEEDABLE_STAGE_DIRS (e.g. a
        leftover scratch directory, or a future stage this list of
        constants hasn't been updated for yet) must not be copied — only
        the explicitly recognised stage names are eligible."""
        src = tmp_path / "old_run"
        dst = tmp_path / "new_run"
        src.mkdir()
        (src / "length_filter").mkdir()
        _touch(src / "length_filter" / "OG0001")
        (src / "some_unrecognised_scratch_dir").mkdir()
        _touch(src / "some_unrecognised_scratch_dir" / "junk")

        _seed_fork_source(str(src), str(dst))

        assert (dst / "length_filter").is_dir()
        assert not (dst / "some_unrecognised_scratch_dir").exists()

    def test_all_known_stage_names_present(self, tmp_path):
        """Sanity check that every stage this test suite exercises above
        is actually a member of the real constant, not a typo'd stand-in
        — keeps this test file honest against the real list used by
        _run_single_list()."""
        for stage in ("length_filter", "group_filter", "alignment",
                     "alignment_trimmed", "ident_alignment", "tables"):
            assert stage in _FORK_SEEDABLE_STAGE_DIRS

    def test_destination_parent_need_not_pre_exist(self, tmp_path):
        """Covers the two-list-mode case: dest_listdir itself
        (e.g. outdir/list1) does not exist yet at seeding time — only the
        top-level outdir does. shutil.copytree must create the full
        destination path itself."""
        src = tmp_path / "old_run" / "list1"
        dst = tmp_path / "new_run" / "list1"   # 'new_run' does not exist yet
        src.mkdir(parents=True)
        (src / "length_filter").mkdir()
        _touch(src / "length_filter" / "OG0001")

        _seed_fork_source(str(src), str(dst))

        assert (dst / "length_filter" / "OG0001").is_file()

    def test_raises_if_fork_source_does_not_exist(self, tmp_path):
        with pytest.raises(EnhydraIOError, match="not found"):
            _seed_fork_source(str(tmp_path / "nonexistent"), str(tmp_path / "new_run"))

    def test_label_included_in_error_message(self, tmp_path):
        with pytest.raises(EnhydraIOError, match="List 1"):
            _seed_fork_source(
                str(tmp_path / "nonexistent"), str(tmp_path / "new_run"),
                label="List 1",
            )

    def test_nested_file_contents_preserved(self, tmp_path):
        """Confirms this is a genuine recursive copy, not just top-level
        filenames — a stage directory's own multi-level contents (e.g.
        alignment_trimmed/ plus its colnumbering sidecar files living in
        a sibling directory) must round-trip byte-for-byte."""
        src = tmp_path / "old_run"
        dst = tmp_path / "new_run"
        src.mkdir()
        colnum_dir = src / "alignment_trimmed_colnumbering"
        colnum_dir.mkdir()
        (colnum_dir / "OG0001.aln.colnumbering").write_text("#ColumnsMap\t0, 1, 2\n")

        _seed_fork_source(str(src), str(dst))

        copied = dst / "alignment_trimmed_colnumbering" / "OG0001.aln.colnumbering"
        assert copied.read_text() == "#ColumnsMap\t0, 1, 2\n"
