import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from organize_10x import organize


class OrganizeTest(unittest.TestCase):
    def test_preview_and_apply_with_spaces(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "GSM GU0744_P matrix.mtx.gz"
            source.write_bytes(b"original")
            ignored = root / "notes.txt"
            ignored.write_text("keep")
            self.assertEqual(len(organize(root)), 1)
            self.assertTrue(source.exists())
            organize(root, apply=True)
            self.assertFalse(source.exists())
            self.assertEqual((root / "GU0744_P" / source.name).read_bytes(), b"original")
            self.assertEqual(ignored.read_text(), "keep")
            self.assertEqual(organize(root), [])

    def test_collision_aborts_before_any_moves(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            earlier = root / "GU0700_P_a.txt"
            later = root / "GU0744_T_b.txt"
            earlier.write_text("first")
            later.write_text("source")
            target_dir = root / "GU0744_T"
            target_dir.mkdir()
            (target_dir / later.name).write_text("existing")
            with self.assertRaises(FileExistsError):
                organize(root, apply=True)
            self.assertTrue(earlier.exists())
            self.assertEqual(later.read_text(), "source")
            self.assertEqual((target_dir / later.name).read_text(), "existing")

    def test_symlinks_and_invalid_sample_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            actual = root / "notes.txt"
            actual.write_text("keep")
            (root / "GU0744_T_link.txt").symlink_to(actual)
            (root / "GU0744_TEST.txt").write_text("not a sample")
            self.assertEqual(organize(root, apply=True), [])
            source = root / "GU0744_P_matrix.mtx"
            source.write_text("matrix")
            (root / "GU0744_P").symlink_to(root, target_is_directory=True)
            with self.assertRaises(ValueError):
                organize(root, apply=True)
            self.assertTrue(source.exists())


if __name__ == "__main__":
    unittest.main()
